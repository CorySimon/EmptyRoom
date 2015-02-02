CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch sm_20 g
NVCCFLAGS = -g -O3 -arch sm_20
LIBS = -lm

all: test

test: test.o Framework.o Forcefield.o
	$(CC) -o $@ $(NVCCLIBS) -L/usr/local/cuda/lib64 -lcurand test.o Framework.o Forcefield.o

Framework.o: Framework.cpp
	$(CC) -c $(NVCCFLAGS) Framework.cpp

Forcefield.o: Forcefield.cpp
	$(CC) -c $(NVCCFLAGS) Forcefield.cpp

clean: 
	rm *.o 