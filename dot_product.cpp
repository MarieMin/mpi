#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <cstdlib>

using namespace std;

#define N 1024

double dot_product(double *a, double *b, int start_idx, int end_idx) {
	int i;
	double product = 0.0;

	for (i = 0; i < end_idx; i++)
	{
		product += a[i] * b[i];
	}
	return product;
}

int main(int argc, char **argv) {

	double *a, *b, *local_a, *local_b;
	double product, local_product;

	double start, stop, seq_dur;
	
	int i;
	int proc_rank, size;
	int tag_a = 1;
	int tag_b = 2;
	int tag_prod = 3;
	
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	const int block_size = N / size;

	local_a = (double*)malloc(block_size * sizeof(double));
	local_b = (double*)malloc(block_size * sizeof(double));

	if (proc_rank == 0)
	{
		printf("\nComm size: %d\n", size);

		a = (double*)malloc(N * sizeof(double));
		b = (double*)malloc(N * sizeof(double));


		for (i = 0; i < N; i++)
		{
			a[i] = rand() % 10;
			b[i] = rand() % 10;
		}

		start = MPI_Wtime();

		// master local copy
		for (int i = 0; i < block_size; i++)
			local_a[i] = a[i];
			local_b[i] = b[i];

		// distributes the portion of array to child processes
		for (i = 1; i < size; i++) {
			MPI_Send(&a[i*block_size],
				block_size,
				MPI_DOUBLE, i, tag_a,
				MPI_COMM_WORLD);

			MPI_Send(&b[i*block_size],
				block_size,
				MPI_DOUBLE, i, tag_b,
				MPI_COMM_WORLD);
		}

		// master portion of scalar product
		product = dot_product(local_a, local_b, 0, block_size);

		// collect local products from other processes 
		double tmp;
		for (i = 1; i < size; i++) {
			MPI_Recv(&tmp, 1, MPI_DOUBLE, i, tag_prod,
				MPI_COMM_WORLD, &status);

			product += tmp;
		}

		stop = MPI_Wtime();
	
		printf("Parallel dot product of %d elements: %f\n", N, product);
		printf("Time: %lf\n", stop - start);
		

		// sequential section
		start = MPI_Wtime();
		double seq_product = dot_product(a, b, 0, N);
		stop = MPI_Wtime();

		printf("Seq dot product of %d elements: %f\n", N, seq_product);
		printf("Time: %f\n\n", stop - start);	

	}
	else {
	
		MPI_Recv(local_a, block_size, MPI_DOUBLE, 0, tag_a,
			MPI_COMM_WORLD, &status);

		MPI_Recv(local_b, block_size, MPI_DOUBLE, 0, tag_b,
			MPI_COMM_WORLD, &status);

		local_product = dot_product(local_a, local_b, 0, block_size);
		
		// sends local result to the root process 
		MPI_Send(&local_product, 1, MPI_DOUBLE,
			0, tag_prod, MPI_COMM_WORLD);
	}


	MPI_Finalize();
	return 0;
}
