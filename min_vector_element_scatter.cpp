#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <limits>

using namespace std;

#define N 32

double payload(double x) {
        return sinh(atan(cbrt(log(exp(pow(tan(asinh(x)),3.0))))));
}

double custom_min(double *a, int start_idx, int end_idx) {
        int i;
        double min = std::numeric_limits<double>::max();

        for (i = start_idx; i < end_idx; i++)
        {
                if (payload(a[i]) < min) { min = a[i]; }
        }
        return min;
}

main(int argc, char* argv[])
{
        double *a, *local_a;
        double min; 
        double local_min;

        double start, stop;

        int i;
        int proc_rank, size;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

        if (proc_rank == 0) {
                a = (double *)calloc(N, sizeof(double));

                for (i = 0; i < N; i++)
                {
                        a[i] = rand() % 100;

                }

                start = MPI_Wtime();
        }

        const int block_size = N / size;

        local_a = (double *)calloc(block_size, sizeof(double));

        MPI_Scatter(a, block_size, MPI_DOUBLE, local_a, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        local_min = custom_min(local_a, 0, block_size);

        free(local_a); 

        MPI_Reduce(&local_min, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

        if (proc_rank == 0) {
                stop = MPI_Wtime();
                printf("Parallel min of %d elements: %f\n", N, min);
                printf("Time: %lf\n", stop - start);

        }

        if (proc_rank == 0) {
                start = MPI_Wtime();
                double seq_min = custom_min(a, 0, N);

                stop = MPI_Wtime();
                printf("Seq dot product of %d elements: %f\n", N, seq_min);
                printf("Time: %f\n\n", stop - start);

                free(a);
        }

        MPI_Finalize();
}

