GCC=gcc

main:
	main.c
	$GCC -o main main.c -lm -ggdb3