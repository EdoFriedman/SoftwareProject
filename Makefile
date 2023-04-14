CC=gcc
CFLAGS=-ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans: spkmeans.o spkmeans.h graph.o graph.h matrix.o matrix.h
	$(CC) -o spkmeans spkmeans.o $(CFLAGS)

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)

matrix.o: matrix.c
	$(CC) -c matrix.c $(CFLAGS)

graph.o: graph.c
	$(CC) -c graph.c $(CFLAGS)