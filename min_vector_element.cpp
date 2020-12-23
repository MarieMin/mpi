#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <cstdlib>
#include <limits>
#include <math.h>

using namespace std;

#define N 1000000000

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

int main(int argc, char **argv) {

        double *a, *local_a;
        double min, local_min;

        double start, stop, seq_dur;

        int i;
        int proc_rank, size;
        int tag_a = 1;
        int tag_min = 2;

        MPI_Status status;
        MPI_Request req;

       MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

        const int block_size = N / size;

        local_a = (double*)malloc(block_size * sizeof(double));

        if (proc_rank == 0)
        {
                printf("\nComm size %d and vecotr size %d \n", size, N);

                a = (double*)malloc(N * sizeof(double));

                for (i = 0; i < N; i++)
                {
                        a[i] = rand() % 100;
                }

                start = MPI_Wtime();

                // master local copy
                for (int i = 0; i < block_size; i++)
                        local_a[i] = a[i];

                // distributes the portion of array to child processes
                for (i = 1; i < size; i++) {
                        MPI_Isend(&a[i*block_size],
                                block_size,
                                MPI_DOUBLE, i, tag_a,
                                MPI_COMM_WORLD, &req);
                        MPI_Wait(&req, &status);
                }
                // master portion of scalar product
                min = custom_min(local_a, 0, block_size);

                // collect local products from other processes 
                double tmp;
                for (i = 1; i < size; i++) {
                        MPI_Irecv(&tmp, 1, MPI_DOUBLE, i, tag_min,
                                MPI_COMM_WORLD, &req);
                        MPI_Wait(&req, &status);

                        if (tmp < min) { min = tmp; } 
                }

                stop = MPI_Wtime();

                printf("Parallel min of %d elements: %f\n", N, min);
                printf("%lf\n", stop - start);


                // sequential section
                start = MPI_Wtime();
                double seq_min = custom_min(a, 0, N);
                stop = MPI_Wtime();

                printf("%f\n", stop - start);
                printf("Seq min: %f\n\n", seq_min);
                free(a);

        }
        else {

                MPI_Irecv(local_a, block_size, MPI_DOUBLE, 0, tag_a,
                        MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);

                local_min = custom_min(local_a, 0, block_size);

                // sends local result to the root process 
                MPI_Isend(&local_min, 1, MPI_DOUBLE,
                        0, tag_min, MPI_COMM_WORLD, &req);
                MPI_Wait(&req, &status);
        }

        free(local_a);
        MPI_Finalize();
        return 0;
}


