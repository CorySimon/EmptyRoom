CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch sm_20 g
NVCCFLAGS = -g -O3 -arch sm_20
LIBS = -lm

all: test

test: test.o Framework.o
	$(CC) -o $@ $(NVCCLIBS) -L/usr/local/cuda/lib64 -lcurand test.o Framework.o 

Framework.o: Framework.cpp
	$(CC) -c $(NVCCFLAGS) Framework.cpp

clean: 
	rm *.o 
