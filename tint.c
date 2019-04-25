#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define A 2.0

double f(double x)
{
	return sqrt(4.0 - x * x);
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		printf("Usage: %s num_intervals.\n", argv[0]);
		return 1;
	}
	MPI_Init(&argc, &argv);
	int rank;
	int i;
	int N = atoi(argv[1]);
	double h = A / N;
	double S = 0.0;
        double S_i = 0.0;
        int P = 0;
        MPI_Comm_size(MPI_COMM_WORLD, &P);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	for (i = rank; i < N; i += P) {
		S_i += h * (f(h * i) + f(h * (i + 1))) / 2.0;
	}
        if (rank == 0) {
                S = S_i;
                for (i = 1; i < P; i++) {
                        MPI_Recv(&S_i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        S += S_i;
                }
                printf("%f\n", S);

        }
        else MPI_Send(&S_i, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
}

