COMPILER = gcc

FLAGS_WARNINGS = -Wextra -Wall -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Waggregate-return -Wformat=2 -Wno-unknown-pragmas 
FLAGS_SETTINGS = -O2 -march=native
FLAGS = $(FLAGS_SETTINGS) $(FLAGS_WARNINGS)
FLAGS_VECTORIZATION = $(FLAGS) -ftree-vectorize
FLAGS_ASSEMBLY = -masm=intel -S -fno-asynchronous-unwind-tables

FILE = main.c
LIB = -lm -lgmp

.PHONY: all
all: run

open:
	nvim main.s

run: clean
	gcc -march=native -O3 test.c sub.c -o test
	gcc -march=native -O3 -S sub.c
	@$(COMPILER) $(FLAGS_VECTORIZATION) -o $@ $(FILE) $(LIB)
	@./$@ $(ALL_ARGUMENTS)

assembly: clean
	$(COMPILER) $(FLAGS_VECTORIZATION) $(FLAGS_ASSEMBLY) -fverbose-asm $(FILE) $(LIB)

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
