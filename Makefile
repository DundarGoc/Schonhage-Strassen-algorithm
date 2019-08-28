COMPILER = gcc

FLAGS_WARNINGS = -Wextra -Wall -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Waggregate-return -Wformat=2 -Wno-unknown-pragmas 
FLAGS = $(FLAGS_WARNINGS) -O2 -march=native -ftree-vectorize
FLAGS_ASSEMBLY = -masm=intel -S -fno-asynchronous-unwind-tables

FILE = main.c
LIB = -lm -lgmp

.PHONY: all
all: run

run: clean
	$(COMPILER) $(FLAGS) $(FLAGS_ASSEMBLY) sub.c
	$(COMPILER) $(FLAGS) $(FLAGS_ASSEMBLY) sub_novec.c

assembly: clean
	$(COMPILER) $(FLAGS) $(FLAGS_ASSEMBLY) -fverbose-asm $(FILE) $(LIB)

valgrind: clean
	reset
	$(COMPILER) $(FLAGS) -march=x86-64 -ggdb3 -Og -o $@ $(FILE) $(LIB)
	valgrind --tool=memcheck --leak-check=full --show-reachable=yes --num-callers=20 --track-fds=yes ./$@

vec: clean
	$(COMPILER) $(FLAGS) -o $@ $(FILE) $(LIB) -fopt-info-vec-all=opt.txt
	#sed -i '/butterfly.c/!d' opt.txt
	vi opt.txt

sanitizer: clean
	clang -fsanitize=address -g $(FILE) $(LIB)
	./a.out

debug: clean
	$(COMPILER) $(FLAGS) -Og -ggdb3 $(FILE) $(LIB)
	gdb a.out

profile: clean
	$(COMPILER) $(FLAGS) -fno-inline -o $@ $(FILE) $(LIB)
	perf record ./$@
	perf report

tidy: clean
	clang-tidy -checks='*, -android-*' $(FILE) -- $(LIB)

cppcheck: clean
	cppcheck --enable=all --suppress=missingIncludeSystem $(FILE)

clean:
	rm -rf gmon.out callgrind* vgcore* opt*.txt data.txt cachegrind.out* perf.data \.#* *~ *.s
	find . -maxdepth 1 -type f -executable -exec rm {} +
	clear
