include <stdio.h>
#include "mpi.h"
#include <stdlib.h>
#include <cstdlib>

using namespace std;

double dot_product(double *a, double *b, int n) {
        int i;
        double product = 0.0;

        for (i = 0; i < n; i++)
        {
                product += a[i] * b[i];
        }
        return product;
}

int main(int argc, char **argv) {

        const int N = 100000;
        int i;

        double a[N];
        double b[N];

        for (i = 0; i < N; i++)
        {
                a[i] = rand() % 10;
                b[i] = rand() % 10;
        }

        double product;
        int proc_rank;
        int size;
        double start, stop, parallel_dur, seq_dur;
        double max_time, min_time, avg_time;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

        if (proc_rank == 0)
        {
                printf("\nComm size: %d\n", size);
        }

        const int block_size = N / size;
        double block_a[block_size];
        double block_b[block_size];

        for (i = 0; i < block_size; i++)
        {
                block_a[i] = a[i + proc_rank * block_size];
                block_b[i] = b[i + proc_rank * block_size];
        }

        //  MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        double local_product;
        local_product = dot_product(block_b, block_a, block_size);
        MPI_Reduce(&local_product, &product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        stop = MPI_Wtime();

        parallel_dur = stop - start;
        /*compute max, min, and average timing statistics*/
        MPI_Reduce(&parallel_dur, &max_time, 1, MPI_DOUBLE,MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&parallel_dur, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
        MPI_Reduce(&parallel_dur, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
        if (proc_rank == 0) {
                printf("Parallel dot product of %d elements: %f\n", N, product);
                avg_time /= size;
                printf("Min: %lf  Max: %lf  Avg:  %lf\n", min_time, max_time,avg_time);
        }

        if (proc_rank == 0) {
                start = MPI_Wtime();
                double seq_product = dot_product(a, b, N);
                stop = MPI_Wtime();

                printf("Seq dot product of %d elements: %f\n", N, seq_product);
                printf("Time: %f\n\n", stop - start);
        }


        MPI_Finalize();
        return 0;
}
