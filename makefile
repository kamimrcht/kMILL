#~ CC=/usr/bin/g++
CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread -fmax-errors=3 -flto
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=kMILL

all: $(EXEC)

kMILL:	main.o	compaction.o	utils.o	readAndSortInputFile.o
	$(CC) -o $@ $^ $(LDFLAGS)

utils.o: 	utils.cpp	utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

readAndSortInputFile.o:	readAndSortInputFile.cpp	readAndSortInputFile.h	compaction.h	utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

compaction.o:	compaction.cpp	compaction.h	utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

main.o:	main.cpp	compaction.h	readAndSortInputFile.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
