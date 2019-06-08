#include <assert.h>
#include <flint/nmod_poly.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include "SSA.h"

#define BENCHMARK_FUNCTION 0
#define NUMBER_OF_ROWS_RESULT 12
#define NUMBER_OF_COLUMNS_RESULT 13
#define LOWER_LIMIT 0.8
#define UPPER_LIMIT 1.2

static double timeElapsed;
static double measurementTimeInSeconds = 0.01;

#define CEIL_DIV(a, b) ((((a) - 1) / (b)) + 1)

/*
   Simple custom timer.
 */
double GetTime(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec * 1e-6;
}

void WriteResultsToFile(double performanceResults[][NUMBER_OF_COLUMNS_RESULT], double modAverage[], double degreeAverage[], double totalAverage, double maxRatio, double minRatio,
                        double maxModRatio, double minModRatio, double maxDegreeRatio, double minDegreeRatio)
{
	FILE *f = fopen("data.txt", "w");
	fprintf(f, "Performance results:\n");
	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		fprintf(f, "%li ", (i + 1) * 5);

		for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
		{
			fprintf(f, "& %.2lf ", performanceResults[i][j]);
		}

		fprintf(f, " \\\\\n");
	}

	fprintf(f, "\nModulo average:\n");
	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		fprintf(f, "%.2lf  ", modAverage[i]);

		if(i != NUMBER_OF_ROWS_RESULT - 1)
		{
			fprintf(f, "&");
		}
	}

	fprintf(f, "\n\nDegree average:\n");
	for(ulong i = 0; i < NUMBER_OF_COLUMNS_RESULT; ++i)
	{
		fprintf(f, "%.2lf  ", degreeAverage[i]);

		if(i != NUMBER_OF_COLUMNS_RESULT - 1)
		{
			fprintf(f, "&");
		}
	}

	fprintf(f, "\n\nTotal Average Ratio: %f", totalAverage);
	fprintf(f, "\nMax Ratio: %f", maxRatio);
	fprintf(f, "\nMin Ratio: %f", minRatio);

	fprintf(f, "\nMax Mod Ratio: %f", maxModRatio);
	fprintf(f, "\nMin Mod Ratio: %f", minModRatio);

	fprintf(f, "\nMax Degree Ratio: %f", maxDegreeRatio);
	fprintf(f, "\nMin Degree Ratio: %f", minDegreeRatio);

	fclose(f);

}

void PrintCharNTimes(char character, ulong numberOfTimes)
{
	for(ulong i = 0; i < numberOfTimes; ++i)
	{
		putchar(character);
	}
}

void PrintInColor(double value)
{
	char red[] = "\x1b[31m";
	char green[] = "\x1b[32m";
	char blue[] = "\x1b[34m";
	char reset[] = "\x1b[0m";

	if(value < LOWER_LIMIT)
	{
		printf("%s%.2lf%s  ", red, value, reset);
	}
	else if(value > UPPER_LIMIT)
	{
		printf("%s%.2lf%s  ", green, value, reset);
	}
	else
	{
		printf("%s%.2lf%s  ", blue, value, reset);
	}
}

void PrintResults(double performanceResult[][NUMBER_OF_COLUMNS_RESULT], double modAverage[], double degreeAverage[], double totalAverage)
{

	printf("\n        ");

	for(ulong i = 0; i < NUMBER_OF_COLUMNS_RESULT; ++i)
	{
		PrintInColor(degreeAverage[i]);
	}

	putchar('\n');
	PrintCharNTimes(' ', 7);
	PrintCharNTimes('_', 6 * NUMBER_OF_COLUMNS_RESULT - 1);
	putchar('\n');

	for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
	{
		PrintInColor(modAverage[i]);
		printf("| ");

		for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
		{
			PrintInColor(performanceResult[i][j]);
		}

		putchar('\n');
	}

	if(BENCHMARK_FUNCTION == 0)
	{
		printf("\nBenchmarking butterfly.");
	}
	else if(BENCHMARK_FUNCTION == 1)
	{
		printf("\nBenchmarking TFT.");
	}
	else
	{
		printf("\nBenchmarking ITFT.");
	}
	printf("\nThe numbers show how much faster the vectorization is.\n");
	printf("Total Average Ratio: ");

	PrintInColor(totalAverage);
	printf("\n\n");

}

void BenchmarkITFT(ulong *poly, ulong polyLength, nmod_t modFLINT)
{
	ulong lgM, m1, m2, m3, pmfLength;

	for(lgM = 1;; ++lgM)
	{
		pmfLength = 1UL << lgM;
		m1 = CEIL_DIV(polyLength, pmfLength / 2);
		m2 = CEIL_DIV(polyLength, pmfLength / 2);
		m3 = m1 + m2 - 1;

		if(m3 <= 2 * pmfLength)
		{
			break;
		}
	}

	ulong lgK = (m3 > pmfLength) ? (lgM + 1) : lgM;
	ulong pmfVectorLength = 1UL << lgK;
	ptrdiff_t pmfLengthWithBias = pmfLength + 1;
	ulong *pmfVector1 = (ulong *)flint_malloc(sizeof(ulong) * pmfLengthWithBias * pmfVectorLength);

	TransformPolynomialToPMFVector(poly, polyLength, 1, pmfLength, modFLINT, pmfVector1);

	double t = GetTime();
	PMFVectorITFT(m3, 0, m3, 0, lgM, lgK, pmfLengthWithBias, pmfVector1, modFLINT);
	timeElapsed += GetTime() - t;

	flint_free(pmfVector1);
}

void BenchmarkTFT(ulong *poly, ulong polyLength, nmod_t modFLINT)
{
	ulong lgM, m1, m2, m3, pmfLength;

	for(lgM = 1;; ++lgM)
	{
		pmfLength = 1UL << lgM;
		m1 = CEIL_DIV(polyLength, pmfLength / 2);
		m2 = CEIL_DIV(polyLength, pmfLength / 2);
		m3 = m1 + m2 - 1;

		if(m3 <= 2 * pmfLength)
		{
			break;
		}
	}

	ulong lgK = (m3 > pmfLength) ? (lgM + 1) : lgM;
	ulong pmfVectorLength = 1UL << lgK;
	ptrdiff_t pmfLengthWithBias = pmfLength + 1;
	ulong *pmfVector1 = (ulong *)flint_malloc(sizeof(ulong) * pmfLengthWithBias * pmfVectorLength);

	TransformPolynomialToPMFVector(poly, polyLength, 1, pmfLength, modFLINT, pmfVector1);

	double t = GetTime();
	PMFVectorTFT(m3, m1, 0, lgK, lgM, pmfLengthWithBias, pmfVector1, modFLINT);
	timeElapsed += GetTime() - t;

	flint_free(pmfVector1);
}

void BenchmarkVectorization(unsigned int modBits, ulong maximumPossibleLength, flint_rand_t state)
{
	assert(modBits >= 2);
	assert(modBits <= 64);
	assert(maximumPossibleLength >= 1);

	ulong mod = n_randbits(state, modBits);

	while(mod % 2 == 0)
	{
		mod = n_randbits(state, modBits);
	}

	ulong polyLength = n_randint(state, maximumPossibleLength) + 1;
	ulong *poly1 = (ulong *)flint_calloc(polyLength, sizeof(ulong));
	ulong *poly2 = (ulong *)flint_calloc(polyLength, sizeof(ulong));
	nmod_t modFLINT;
	nmod_init(&modFLINT, mod);

	for(ulong i = 0; i < polyLength; ++i)
	{
		poly1[i] = n_randlimb(state);
		poly2[i] = n_randlimb(state);
	}

	_nmod_vec_reduce(poly1, poly1, polyLength, modFLINT);
	_nmod_vec_reduce(poly2, poly2, polyLength, modFLINT);

	if(BENCHMARK_FUNCTION == 0)
	{
		double t = GetTime();
		ButterflyInPlace(poly1, poly2, polyLength, modFLINT);
		timeElapsed += GetTime() - t;
	}
	else if(BENCHMARK_FUNCTION == 1)
	{
		BenchmarkTFT(poly1, polyLength, modFLINT);
	}
	else
	{
		BenchmarkITFT(poly1, polyLength, modFLINT);
	}

	flint_free(poly1);
	flint_free(poly2);
}

int main(int argc, char *argv[])
{
	char *ptr;
	if(argc != 2)
	{
		printf("Wrong number of arguments. Setting measurement time to 0.01.\n");
	}
	else
	{
		measurementTimeInSeconds = strtod(argv[1], &ptr);
	}

	FILE *f;
	flint_rand_t state;
	flint_randinit(state);
	flint_randseed(state, time(0), time(0));
	setbuf(stdout, NULL);

	double data[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
	double dataRunOne[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
	ulong numberOfIterationsTable[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT] = {0};
	int currentlyOnRunOne = (f = fopen("data.txt", "r")) == NULL;
	int currentlyOnRunTwo = !currentlyOnRunOne;
	if(currentlyOnRunTwo)
	{
		if(!fread(dataRunOne, sizeof(double), NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT, f))
		{
			printf("Failed to read dataRunOne.\n");
		}

		fclose(f);

		if((f = fopen("numberOfIterations.txt", "r")) != NULL)
		{
			if(!fread(numberOfIterationsTable, sizeof(ulong), NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT, f))
			{
				printf("Failed to read dataRunOne.\n");
			}

			fclose(f);

		}
		else
		{
			printf("\nFailed to open numberOfIterations.txt\n");
		}
	}

	for(unsigned int modBits = 5; modBits <= 60; modBits += 5)
	{
		for(ulong idx = 0; idx < NUMBER_OF_COLUMNS_RESULT; ++idx)
		{
			ulong maximumPossibleLength = 1024 * (1UL << idx / 2);
			if(idx & 1UL)
			{
				maximumPossibleLength += maximumPossibleLength / 2;
			}

			timeElapsed = 0;
			if(currentlyOnRunOne)
			{
				double startTime = GetTime();
				ulong numberOfIterations = 0;
				while(GetTime() - startTime < measurementTimeInSeconds)
				{
					BenchmarkVectorization(modBits, maximumPossibleLength, state);
					++numberOfIterations;
				}
				numberOfIterationsTable[modBits / 5 - 1][idx] = numberOfIterations;
				data[modBits / 5 - 1][idx] = timeElapsed;
			}
			else
			{
				for(ulong i = 0; i < numberOfIterationsTable[modBits / 5 - 1][idx]; ++i)
				{
					BenchmarkVectorization(modBits, maximumPossibleLength, state);
				}
				data[modBits / 5 - 1][idx] = timeElapsed;
			}
		}
	}

	if(currentlyOnRunOne)
	{
		f = fopen("data.txt", "w");
		fwrite(data, sizeof(double), NUMBER_OF_COLUMNS_RESULT * NUMBER_OF_ROWS_RESULT, f);
		fclose(f);

		f = fopen("numberOfIterations.txt", "w");
		fwrite(numberOfIterationsTable, sizeof(ulong), NUMBER_OF_COLUMNS_RESULT * NUMBER_OF_ROWS_RESULT, f);
		fclose(f);
	}
	else
	{
		double performanceResult[NUMBER_OF_ROWS_RESULT][NUMBER_OF_COLUMNS_RESULT];
		double totalAverage = 0;
		double maxRatio = 0;
		double minRatio = INFINITY;
		double modAverage[NUMBER_OF_ROWS_RESULT] = {0};
		double maxModRatio = 0;
		double minModRatio = INFINITY;
		double degreeAverage[NUMBER_OF_COLUMNS_RESULT] = {0};
		double maxDegreeRatio = 0;
		double minDegreeRatio = INFINITY;

		for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
		{
			for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
			{
				performanceResult[i][j] = dataRunOne[i][j] / data[i][j];
				maxRatio = FLINT_MAX(maxRatio, performanceResult[i][j]);
				minRatio = FLINT_MIN(minRatio, performanceResult[i][j]);
			}
		}

		for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
		{
			for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
			{
				modAverage[i] += performanceResult[i][j] / NUMBER_OF_COLUMNS_RESULT;
				maxModRatio = FLINT_MAX(maxModRatio, modAverage[i]);
				minModRatio = FLINT_MIN(minModRatio, modAverage[i]);
			}
		}

		for(ulong i = 0; i < NUMBER_OF_COLUMNS_RESULT; ++i)
		{
			for(ulong j = 0; j < NUMBER_OF_ROWS_RESULT; ++j)
			{
				degreeAverage[i] += performanceResult[j][i] / NUMBER_OF_ROWS_RESULT;
				maxDegreeRatio = FLINT_MAX(maxDegreeRatio, degreeAverage[i]);
				minDegreeRatio = FLINT_MIN(minDegreeRatio, degreeAverage[i]);
			}
		}

		for(ulong i = 0; i < NUMBER_OF_ROWS_RESULT; ++i)
		{
			for(ulong j = 0; j < NUMBER_OF_COLUMNS_RESULT; ++j)
			{
				totalAverage += performanceResult[i][j] / (double)(NUMBER_OF_ROWS_RESULT * NUMBER_OF_COLUMNS_RESULT);
			}
		}

		PrintResults(performanceResult, modAverage, degreeAverage, totalAverage);
		WriteResultsToFile(performanceResult, modAverage, degreeAverage, totalAverage, maxRatio, minRatio, maxModRatio, minModRatio, maxDegreeRatio, minDegreeRatio);
	}

	flint_randclear(state);
}
