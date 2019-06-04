MEASUREMENT_TIME_IN_SECONDS = 0.0001
CHOOSE_FUNCTION_1 = 1  #0 - SSA, 1 - FLINT, 2 - zn_poly
CHOOSE_FUNCTION_2 = 0
ALL_ARGUMENTS = $(MEASUREMENT_TIME_IN_SECONDS) $(CHOOSE_FUNCTION_1) $(CHOOSE_FUNCTION_2)

COMPILER = gcc

FLAGS_WARNINGS = -Wextra -Wall -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Waggregate-return -Wformat=2 -Wno-unknown-pragmas 
FLAGS_SETTINGS = -O2 -march=native -pipe -DNDEBUG
FLAGS = $(FLAGS_SETTINGS) $(FLAGS_WARNINGS)
FLAGS_VECTORIZATION = $(FLAGS) -ftree-vectorize
FLAGS_ASSEMBLY = -masm=intel -S -fno-asynchronous-unwind-tables

#FILE = main.c tuning.c SSA.c
FILE = benchmarkVectorization.c tuning.c SSA.c

ZN_POLY = -I/home/dundar/local_zn_poly/include -L/home/dundar/local_zn_poly/lib -lzn_poly
FLINT = -Wl,-rpath=/home/dundar/local_flint/lib -I/home/dundar/local_flint/include -L/home/dundar/local_flint/lib -lflint
LIB = $(FLINT) $(ZN_POLY) -lm -lgmp


.PHONY: all
all: benchmarkVectorization

benchmarkMultiplication: clean
	@$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB)
	@./$@ $(ALL_ARGUMENTS)

benchmarkVectorization: clean
	@$(COMPILER) $(FLAGS) -o $@ $(FILE) $(LIB)
	./$@ $(MEASUREMENT_TIME_IN_SECONDS)
	@$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB)
	./$@ $(MEASUREMENT_TIME_IN_SECONDS)

assembly: clean
	$(COMPILER) $(FLAGS_VECTORIZATION) $(FLAGS_ASSEMBLY) -fverbose-asm	 $(FILE) $(LIB)

valgrind: clean
	reset
	$(COMPILER) $(FLAGS_VECTORIZATION) -march=x86-64 -ggdb3 -Og -o $@ $(FILE) $(LIB)
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --num-callers=20 --track-fds=yes ./$@ $(ALL_ARGUMENTS)

vec: clean
	$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB) -fopt-info-vec-all=opt.txt
	#sed -i '/butterfly.c/!d' opt.txt
	vi opt.txt

sanitizer: clean
	clang -fsanitize=address -g $(FILE) $(LIB)
	./a.out $(ALL_ARGUMENTS)

debug: clean
	$(COMPILER) $(FLAGS_VECTORIZATION) -Og -ggdb3 $(FILE) $(LIB)
	gdb a.out

profile: clean
	$(COMPILER) $(FLAGS_VECTORIZATION) -fno-inline -o $@ $(FILE) $(LIB)
	perf record ./$@ $(ALL_ARGUMENTS)
	perf report

tidy: clean
	clang-tidy -checks='*, -android-*' $(FILE) -- $(LIB)

cppcheck: clean
	cppcheck --enable=all --suppress=missingIncludeSystem $(FILE)

clean:
	rm -rf gmon.out callgrind* vgcore* opt*.txt data.txt cachegrind.out* perf.data numberOfIterations.txt \.#* *~ *.s
	find . -maxdepth 1 -type f -executable -exec rm {} +
	clear