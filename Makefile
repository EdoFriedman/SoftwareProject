CC=gcc
CFLAGS=-ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans: spkmeans.o spkmeans.h matrix.o matrix.h
	$(CC) -o spkmeans spkmeans.o matrix.o $(CFLAGS)

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)

matrix.o: matrix.c
	$(CC) -c matrix.c $(CFLAGS)
