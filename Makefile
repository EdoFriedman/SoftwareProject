CC=gcc
CFLAGS=-ansi -Wall -Wextra -Werror -pedantic-errors -lm

spkmeans.o: spkmeans.c
	$(CC) -c spkmeans.c $(CFLAGS)

spkmeans: spkmeans.o spkmeans.h
	$(CC) -o spkmeans spkmeans.o $(CFLAGS)