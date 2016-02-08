CC=/usr/bin/g++
#CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread
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

kMILL:   main.o  
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp  
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
