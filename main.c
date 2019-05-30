#include <sys/time.h>
#include <time.h>
#include <zn_poly/zn_poly.h>
#include "include/SSA.h"

/*
   Choose benchmark style
 */
#define BENCHMARK_SHOUP_STYLE 1
#define BENCHMARK_LOWER_THRESHOLD 0.8
#define BENCHMARK_UPPER_THRESHOLD 1.2

/*
   These settings are only relevant if BENCHMARK_SHOUP_STYLE = 0.
 */
#define LENGTH_BITS_MINIMUM 1
#define LENGTH_BITS_MAXIMUM 10
#define MODULO_BITS_MINIMUM 2
#define MODULO_BITS_MAXIMUM 64
#define MODULO_BITS_STEP_SIZE 5

/*
   Simple custom timer.
 */
double GetTime(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec * 1e-6;
}

/*
   Global variables.
 */
static double timer1 = 0;
static double timer2 = 0;
static double measurementTimeInSeconds = 0.01;
static int function1 = 1;
static int function2 = 0;

const char *ReturnFunctionName(ulong functionNumber)
{
	switch(functionNumber)
	{
	 case 0:
		 return "SSA";
	 case 1:
		 return "FLINT";
	 default:
		 return "zn_poly";
	}
}

void WriteResultsToFile(ulong numberOfRows, ulong numberOfColumns, double **data)
{
	FILE *f = fopen("data.txt", "w");
	for(ulong i = 0; i < numberOfRows; ++i)
	{
		fprintf(f, "%li ", (i + 1) * 5);

		for(ulong j = 0; j < numberOfColumns; ++j)
		{
			fprintf(f, "& %.2lf ", data[i][j]);

			if(j == numberOfColumns - 1)
			{
				fprintf(f, " \\\\");
			}
		}

		if(i != numberOfRows - 1)
		{
			fprintf(f, "\n");
		}
	}

	fclose(f);
}

void PrintInColor(double value)
{
	char red[] = "\x1b[31m";
	char green[] = "\x1b[32m";
	char blue[] = "\x1b[34m";
	char reset[] = "\x1b[0m";

	if(value < BENCHMARK_LOWER_THRESHOLD)
	{
		printf("%s%.2lf%s", red, value, reset);
	}
	else if(value > BENCHMARK_UPPER_THRESHOLD)
	{
		printf("%s%.2lf%s", green, value, reset);
	}
	else
	{
		printf("%s%.2lf%s", blue, value, reset);
	}
}

void PrintResults(int numberOfRows, int numberOfColumns, double **data, ulong *maximumPossibleLengthAllValues)
{
	double totalAverage = 0;
	for(int i = 0; i < numberOfRows; ++i)
	{
		for(int j = 0; j < numberOfColumns; ++j)
		{
			totalAverage += data[i][j] / (double)(numberOfRows * numberOfColumns);
		}
	}

	char messageMod1[] = "modulo";
	char messageMod2[] = "#bits";
	char messageLength[] = "Max. polynomial length";

	int WIDTH_OF_DATA_POINT = 4;
	int WIDTH_OF_MODULO_BITS = 2;

	int leftMarginAll = 2;
	int spaceAfterModuloText = 2;
	int widthOfModText = FLINT_MAX(sizeof(messageMod1), sizeof(messageMod2));
	int spaceAfterVerticalBar = 1;
	int spaceBeforeVerticalBar = 1;
	int spaceBetweenData = 3;

	int leftMarginHorizontalBar = leftMarginAll + WIDTH_OF_MODULO_BITS + spaceBeforeVerticalBar + widthOfModText + spaceAfterModuloText;
	int horizontalBarLength = WIDTH_OF_DATA_POINT * numberOfColumns + spaceAfterVerticalBar + spaceBetweenData * (numberOfColumns - 1);
	int leftMarginMessageLength = (leftMarginHorizontalBar + horizontalBarLength) / 2;
	int leftMarginMaxPolyLength = leftMarginHorizontalBar + 1 + spaceAfterVerticalBar;

	printf("\n%*s%s", leftMarginMessageLength, "", messageLength);

	printf("\n\n%*s", leftMarginMaxPolyLength, "");
	for(int i = 0; i < numberOfColumns; ++i)
	{
		printf("%lu%*s", maximumPossibleLengthAllValues[i], spaceBetweenData + WIDTH_OF_DATA_POINT - n_sizeinbase(maximumPossibleLengthAllValues[i], 10), "");
	}

	printf("\n%*s+", leftMarginHorizontalBar, "");
	for(int i = 0; i < horizontalBarLength; ++i)
	{
		printf("-");
	}

	ulong moduloBitsMinimum = BENCHMARK_SHOUP_STYLE ? 5 : MODULO_BITS_MINIMUM;
	ulong moduloBitsStepSize = BENCHMARK_SHOUP_STYLE ? 5 : MODULO_BITS_STEP_SIZE;
	for(int i = 0, moduloBits = moduloBitsMinimum; i < numberOfRows; ++i, moduloBits += moduloBitsStepSize)
	{
		printf("\n%*s", leftMarginAll, "");

		if(i == numberOfRows / 2 - 1)
		{
			printf("%s%*s", messageMod1, widthOfModText - (int)sizeof(messageMod1) + 1 + spaceAfterModuloText, "");
		}
		else if(i == numberOfRows / 2)
		{
			printf("%s%*s", messageMod2, widthOfModText - (int)sizeof(messageMod2) + 1 + spaceAfterModuloText, "");
		}
		else
		{
			printf("%*s", widthOfModText + spaceAfterModuloText, "");
		}

		printf("%2i", moduloBits);
		printf("%*s|%*s", spaceBeforeVerticalBar, "", spaceAfterVerticalBar, "");

		for(int j = 0; j < numberOfColumns; ++j)
		{
			PrintInColor(data[i][j]);
			printf("%*s", spaceBetweenData, "");
		}
	}

	printf("\n\nThe numbers show how much faster %s is compared to %s. ", ReturnFunctionName(function2), ReturnFunctionName(function1));
	printf("Measurement time set to %.2f seconds.\nAverage performance increase: ", measurementTimeInSeconds);
	PrintInColor(totalAverage);

	printf("\n\n");
}

void BenchmarkMultiplication(ulong mod, ulong poly1Length, ulong poly2Length, flint_rand_t state)
{
	ulong polyResultLength = poly1Length + poly2Length - 1;
	ulong *poly1 = (ulong *)flint_calloc(poly1Length, sizeof(ulong));
	ulong *poly2 = (ulong *)flint_calloc(poly2Length, sizeof(ulong));
	ulong *polyResult1 = (ulong *)flint_calloc(polyResultLength, sizeof(ulong));
	ulong *polyResult2 = (ulong *)flint_calloc(polyResultLength, sizeof(ulong));
	nmod_t modFLINT;
	nmod_init(&modFLINT, mod);
	zn_mod_t modZN;
	zn_mod_init(modZN, mod);

	for(ulong i = 0; i < poly1Length; ++i)
	{
		poly1[i] = n_randlimb(state);
	}

	for(ulong i = 0; i < poly2Length; ++i)
	{
		poly2[i] = n_randlimb(state);
	}

	_nmod_vec_reduce(poly1, poly1, poly1Length, modFLINT);
	_nmod_vec_reduce(poly2, poly2, poly2Length, modFLINT);

	double elapsedTime = GetTime();
	if(function1 == 0)
	{
		SSA(polyResult1, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(function1 == 1)
	{
		_nmod_poly_mul(polyResult1, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(function1 == 2)
	{
		zn_array_mul(polyResult1, poly1, poly1Length, poly2, poly2Length, modZN);
	}
	elapsedTime = GetTime() - elapsedTime;
	timer1 += elapsedTime;

	elapsedTime = GetTime();
	if(function2 == 0)
	{
		SSA(polyResult2, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(function2 == 1)
	{
		_nmod_poly_mul(polyResult2, poly1, poly1Length, poly2, poly2Length, modFLINT);
	}
	else if(function2 == 2)
	{
		zn_array_mul(polyResult2, poly1, poly1Length, poly2, poly2Length, modZN);
	}
	elapsedTime = GetTime() - elapsedTime;
	timer2 += elapsedTime;


	assert(_nmod_vec_equal(polyResult1, polyResult2, polyResultLength));

	flint_free(poly1);
	flint_free(poly2);
	flint_free(polyResult1);
	flint_free(polyResult2);
}

int main(int argc, char *argv[])
{
	assert(MODULO_BITS_MINIMUM >= 2);
	assert(MODULO_BITS_MAXIMUM <= 64);
	assert(MODULO_BITS_MINIMUM <= MODULO_BITS_MAXIMUM);
	assert(LENGTH_BITS_MINIMUM <= LENGTH_BITS_MAXIMUM);
	assert(LENGTH_BITS_MINIMUM >= 1);

	if(argc != 4)
	{
		printf("Wrong number of arguments. Setting measurement time to 0.01, function 1 to SSA and function 2 to FLINT.\n");
	}
	else
	{
		measurementTimeInSeconds = atof(argv[1]);
		function1 = atoi(argv[2]);
		function2 = atoi(argv[3]);
	}

	double timer = GetTime();
	setbuf(stdout, NULL);

	flint_rand_t state;
	flint_randinit(state);
	flint_randseed(state, time(0), time(0));

	ulong numberOfRows = BENCHMARK_SHOUP_STYLE ? 12 : (MODULO_BITS_MAXIMUM - MODULO_BITS_MINIMUM) / MODULO_BITS_STEP_SIZE + 1;
	ulong numberOfColumns = BENCHMARK_SHOUP_STYLE ? 13 : LENGTH_BITS_MAXIMUM - LENGTH_BITS_MINIMUM + 1;

	double *performanceResultsEntries = (double *) malloc(sizeof(double) * numberOfRows * numberOfColumns);
	double **performanceResults = (double **) malloc(sizeof(double *) * numberOfRows);
	for(ulong ix = 0, jx = 0; ix < numberOfRows; ++ix, jx += numberOfColumns)
	{
		performanceResults[ix] = performanceResultsEntries + jx;
	}

	ulong moduloBitsMinimum = BENCHMARK_SHOUP_STYLE ? 5 : MODULO_BITS_MINIMUM;
	ulong moduloBitsMaximum = BENCHMARK_SHOUP_STYLE ? 60 : MODULO_BITS_MAXIMUM;
	ulong moduloBitsStepSize = BENCHMARK_SHOUP_STYLE ? 5 : MODULO_BITS_STEP_SIZE;
	ulong *maximumPossibleLengthAllValues = (ulong *)flint_calloc(numberOfColumns, sizeof(ulong));

	if(BENCHMARK_SHOUP_STYLE)
	{
		for(ulong column = 0; column < numberOfColumns; ++column)
		{
			long maximumPossibleLength = 1024L << column / 2;

			if(column & 1)
			{
				maximumPossibleLength += maximumPossibleLength / 2;
			}

			maximumPossibleLengthAllValues[column] = maximumPossibleLength;
		}
	}
	else
	{
		for(ulong lengthBits = LENGTH_BITS_MINIMUM, column = 0; lengthBits <= LENGTH_BITS_MAXIMUM; ++lengthBits, ++column)
		{
			maximumPossibleLengthAllValues[column] = (1UL << lengthBits) - 1;
		}
	}



	for(ulong modBits = moduloBitsMinimum, row = 0; modBits <= moduloBitsMaximum; modBits += moduloBitsStepSize, ++row)
	{
		for(ulong column = 0; column < numberOfColumns; ++column)
		{
			long maximumPossibleLength = maximumPossibleLengthAllValues[column];
			double timeElapsed = 0;
			timer1 = 0;
			timer2 = 0;
			while(timeElapsed <= measurementTimeInSeconds)
			{
				ulong poly1Length = n_randint(state, maximumPossibleLength) + 1;
				ulong poly2Length = n_randint(state, maximumPossibleLength) + 1;
				ulong mod = n_randbits(state, modBits);

				while(mod % 2 == 0)
				{
					mod = n_randbits(state, modBits);
				}

				ulong polyBigLength = FLINT_MAX(poly1Length, poly2Length);
				ulong polySmallLength = FLINT_MIN(poly1Length, poly2Length);

				double t = GetTime();
				BenchmarkMultiplication(mod, polyBigLength, polySmallLength, state);
				timeElapsed += GetTime() - t;
			}
			performanceResults[row][column] = timer1 / timer2;
		}
	}

	WriteResultsToFile(numberOfRows, numberOfColumns, performanceResults);
	double totalProgramTime = GetTime() - timer;
	printf("Total program time: %f\n", totalProgramTime);
	PrintResults(numberOfRows, numberOfColumns, performanceResults, maximumPossibleLengthAllValues);


	flint_free(performanceResults);
	flint_free(performanceResultsEntries);
	flint_free(maximumPossibleLengthAllValues);
	flint_randclear(state);
}
