CC = nvcc
MPCC = nvcc
OPENMP = 
CFLAGS = -O3 -arch sm_20 g
NVCCFLAGS = -g -O3 -arch sm_20
LIBS = -lm

all: writegrid henry

test: test.o Framework.o Forcefield.o
	$(CC) -o $@ $(NVCCLIBS) -L/usr/local/cuda/lib64 -lcurand test.o Framework.o Forcefield.o

writegrid: writegrid.o Framework.o Forcefield.o
	$(CC) -o $@ $(NVCCLIBS) -L/usr/local/cuda/lib64 writegrid.o Framework.o Forcefield.o

henry: henry.o Framework.o Forcefield.o
	$(CC) -o $@ $(NVCCLIBS) -L/usr/local/cuda/lib64 henry.o Framework.o Forcefield.o

Framework.o: Framework.cc
	$(CC) -c $(NVCCFLAGS) Framework.cc

Forcefield.o: Forcefield.cc
	$(CC) -c $(NVCCFLAGS) Forcefield.cc

writegrid.o: writegrid.cu 
	$(CC) -c $(NVCCFLAGS) writegrid.cu

henry.o: henry.cu 
	$(CC) -c $(NVCCFLAGS) henry.cu

clean: 
	rm *.o 
