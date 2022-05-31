
CFLAGS = `pkg-config --cflags gsl` -O3
LFLAGS = `pkg-config --libs gsl`

all: seq para

para:
	mpicc $(CFLAGS) para.c $(LFLAGS) -o para

seq:
	gcc $(CFLAGS) seq.c $(LFLAGS) -o seq
	
clean:
	rm para seq
