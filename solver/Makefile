CC = g++
CFLAGS = -Wall
objects = main.o solver.o matrix_operations.o flux_BC.o quad.o residuals.o readgrid.o

all: output
output: $(objects) DG.h  
	$(CC) $(CFLAG) -o main.exe $(objects) -O2 -larmadillo
clean:
	rm *.o main.exe
