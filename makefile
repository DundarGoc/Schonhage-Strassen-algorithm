all:
	gcc -march=native -O3 test.c sub.c -o test
	gcc -march=native -O3 -S sub.c
