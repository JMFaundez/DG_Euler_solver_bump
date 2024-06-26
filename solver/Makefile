CC = g++
CFLAGS = -Wall
objects = main.o solver.o matrix_operations.o flux_BC.o quad.o residuals.o readgrid.o
all: output

output: $(objects) DG.h  
	$(CC) $(CFLAGS) -o main.exe $(objects) -O2 -larmadillo -fopenmp
main.o: main.cpp DG.h 
	$(CC) $(CFLAGS) -c main.cpp -fopenmp -larmadillo
residuals.o: residuals.cpp DG.h 
	$(CC) $(CFLAGS) -c residuals.cpp -fopenmp -larmadillo
clean:
	rm *.o main.exe
