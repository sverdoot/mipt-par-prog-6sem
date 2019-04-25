#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define L 1.0
#define C 1.0
#define T1 1.0
#define T2 2.0

int main(int argc, char **argv)
{
	if (argc < 3) {
		printf("Usage: %s N T\n", argv[0]);
	}
	int N = atoi(argv[1]);
	double T = atof(argv[2]);
	double h = L / (N - 1);
	double dt = 0.3 * h * h / C;
	int steps = T / dt;
	int i, j;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	/* Определение длины участка каждого процесса */
        int n = N / size;
	/* N = x*n + y*(n+1),  x+y=size
	   -> y = N - size*n	*/
	int y = N - size*n;
	int x = size - y;
	int len = 0;
	if (rank < x) {
        	len = n + 2;
	}
	else {
		len = n + 3;
	}
	if (rank == 0 || rank == size - 1)
                len -= 1;
        if (size == 1)
                len = N;
        /* Выделение памяти. */
        double *u1 = (double*)calloc(n + 3, sizeof(double));
        double *u2 = (double*)calloc(n + 3, sizeof(double));
        double *t = u1;

        /* Граничные условия. */
        if (rank == 0)
                u2[0] = u1[0] = T1;
        if (rank == size - 1)
                u2[len - 1] = u1[len - 1] = T2;


        /* Цикл интегрирования. */
        for (i = 0; i < steps; i++) {
                /* Начало обмена. */
                 if (size != 1) {
                        if (rank % 2 == 0) {
                               if (rank != 0) {
                                        MPI_Send(&u1[1], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                                        MPI_Recv(&u1[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                }
                                if (rank != size - 1) {
                                        MPI_Recv(&u1[len-1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                        MPI_Send(&u1[len-2], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                                }
                        }
                        else {
                                if (rank != size - 1) {
                                        MPI_Recv(&u1[len-1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                        MPI_Send(&u1[len-2], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                                }
                                if (rank != 0) {
                                        MPI_Send(&u1[1], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                                        MPI_Recv(&u1[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                                }
			}
		}
                /* Конец обмена. */
                for (j = 1; j < len - 1; j++) {
                        u2[j] = u1[j] + 0.3 * (u1[j-1] - 2 * u1[j] + u1[j+1]);
                }
                t = u1;
                u1 = u2;
        u2 = t;
        }

        if (rank == 0) {
		int s = 0;
                for (j = 0; j < len - 1; j++) {
			printf("%f %f\n", h * s, u1[j]);
			s += 1;
		}
                for (i = 1; i < size; i++) {
			int sz = (i < x) ? n + 2 : n + 3;
			if (i == size - 1)
				sz -= 1;
                        MPI_Recv(u2, sz, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for (j = 1; j < sz - 1; j++) {
				printf("%f %f\n", h * s, u2[j]);
				s += 1;
                        }
                }
		printf("%f %f\n", h * s, T2);
        }
        else {
		MPI_Send(u1, len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
        free(u2);
	free(u1);
        MPI_Finalize();
        return 0;
}


