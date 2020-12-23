#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 1000

int main(int argc, char *argv[]) {

        double start_p, stop_p, start_seq, stop_seq;
        int proc_rank, size;
        int tag = 0;

        char* message = (char *) malloc(N * sizeof(char));
        for (int i = 0; i < N; i++)
                message[i] = 'i';

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        MPI_Status status;

       if (proc_rank < 2) {

         if (proc_rank == 0){ 
                start_p = MPI_Wtime();
                printf("N = %d\n", N);
                MPI_Send(message, N, MPI_CHAR, 1, 0, MPI_COMM_WORLD);
          }
        for (int j = 0; j < 100; j++){

                int src_proc = (proc_rank + 1) %2;
                MPI_Recv(message, N, MPI_CHAR, src_proc, 0, MPI_COMM_WORLD, &status);
        }

        if (proc_rank == 0){ 
                stop_p = MPI_Wtime();
                printf("%f\n", stop_p - start_p);
        }

        }
        MPI_Finalize();
        free(message);
        return 0;
}
