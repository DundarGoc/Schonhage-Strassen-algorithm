MEASUREMENT_TIME_IN_SECONDS = 0.1
CHOOSE_FUNCTION_1 = 1  #0 - SSA, 1 - FLINT, 2 - zn_poly
CHOOSE_FUNCTION_2 = 0
ALL_ARGUMENTS = $(MEASUREMENT_TIME_IN_SECONDS) $(CHOOSE_FUNCTION_1) $(CHOOSE_FUNCTION_2)

COMPILER = gcc

FLAGS_WARNINGS = -Wextra -Wall -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Waggregate-return -Wformat=2 -Wno-unknown-pragmas 
FLAGS_SETTINGS = -O2 -march=native -pipe 
FLAGS = $(FLAGS_SETTINGS) $(FLAGS_WARNINGS)
FLAGS_VECTORIZATION = $(FLAGS) -ftree-vectorize
FLAGS_ASSEMBLY = -masm=intel -S -fno-asynchronous-unwind-tables


#FILE = mainOneFunction.c
#FILE_1 = $(FILE)

FILE = main.c tuning.c SSA.c

#FILE = benchmarkButterfly.c butterfly.c
#FILE_1 = benchmarkButterfly.c butterflyOriginal.c

ZN_POLY = -I/home/dundar/local_zn_poly/include -L/home/dundar/local_zn_poly/lib -lzn_poly
FLINT = -Wl,-rpath=/home/dundar/local_flint/lib -I/home/dundar/local_flint/include -L/home/dundar/local_flint/lib -lflint
LIB = $(FLINT) $(ZN_POLY) -lm -lgmp


.PHONY: all
all: run

run: initialSetup
	@$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB)
	@./$@ $(ALL_ARGUMENTS)

benchmark: initialSetup
	$(COMPILER) $(FLAGS) -o $@ $(FILE_1) $(LIB)
	./$@ $(MEASUREMENT_TIME_IN_SECONDS) $(CHOOSE_FUNCTION_1)
	$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB)
	./$@ $(MEASUREMENT_TIME_IN_SECONDS) $(CHOOSE_FUNCTION_2)

assembly: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) $(FLAGS_ASSEMBLY) -fverbose-asm	 $(FILE) $(LIB)

assemblyCompare: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) $(FLAGS_ASSEMBLY) $(FILE_1) $(LIB)
	$(COMPILER) $(FLAGS_VECTORIZATION) $(FLAGS_ASSEMBLY) $(FILE) $(LIB)
	sed -i '/\.[a-z]/d' butterfly.s butterflyOriginal.s
	clear
	diff -y -W 120 butterflyOriginal.s butterfly.s

valgrind: initialSetup
	reset
	$(COMPILER) $(FLAGS_VECTORIZATION) -march=x86-64 -ggdb3 -Og -o $@ $(FILE) $(LIB)
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --num-callers=20 --track-fds=yes ./$@ $(ALL_ARGUMENTS)

vec: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB) -fopt-info-vec-all=opt.txt
	#sed -i '/butterfly.c/!d' opt.txt
	vi opt.txt

vecCompare: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE_1) $(LIB) -fopt-info-vec-missed=opt1.txt
	$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB) -fopt-info-vec-missed=opt2.txt
	sed -i '/butterflyOriginal.c/!d' opt1.txt
	sed -e s/butterflyOriginal.c://g -i opt1.txt
	sed -i '/butterfly.c/!d' opt2.txt
	sed -e s/butterfly.c://g -i opt2.txt
	reset
	diff -y -W 150 opt1.txt opt2.txt

sanitizer: initialSetup
	clang -fsanitize=address -g $(FILE) $(LIB)
	./a.out $(ALL_ARGUMENTS)

debug: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) -Og -ggdb3 $(FILE) $(LIB)
	gdb a.out

profile: initialSetup
	$(COMPILER) $(FLAGS_VECTORIZATION) -fno-inline -o $@ $(FILE) $(LIB)
	perf record ./$@ $(ALL_ARGUMENTS)
	perf report

tidy: initialSetup
	clang-tidy -checks='*, -android-*' $(FILE) -- $(LIB)

cppcheck:
	cppcheck --enable=all --suppress=missingIncludeSystem $(FILE)

convert: initialSetup
ifneq ("$(wildcard *.c)","")
	for f in *.c; do mv "$$f" $$(echo "$$f" | sed 's/\.c/\.cpp/g'); done
	@printf "\n\nSuccessfully changed .c files to .cpp \n\n"
else
	for f in *.cpp; do mv "$$f" $$(echo "$$f" | sed 's/\.cpp/\.c/g'); done
	@printf "\n\nSuccessfully changed .cpp files to .c \n\n"
endif

initialSetup: check changeNames clean

check:
ifneq ("$(wildcard *.c)","")
ifneq ("$(wildcard *.cpp)","")
	$(error Both .c and .cpp files in directory. Please choose one of them)
endif
endif

changeNames:
ifeq ("$(wildcard *.c)","")
	$(eval FILE = $(FILE:.c=.cpp))
	$(eval FILE_1 = $(FILE_1:.c=.cpp))
	$(eval LIB += -lstdc++)
endif

clean:
	rm -rf gmon.out callgrind* vgcore* opt*.txt data.txt cachegrind.out* perf.data numberOfIterations.txt \.#* *~ *.s
	find . -maxdepth 1 -type f -executable -exec rm {} +
	rm -rf .vscode/ipch
	clear