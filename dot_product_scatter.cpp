#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <cstdlib>

using namespace std;

#define N 1024

double dot_product(double *a, double *b, int start_idx, int end_idx) {
	int i;
	double product = 0.0;

	for (i = start_idx; i < end_idx; i++)
	{
		product += a[i] * b[i];
	}
	return product;
}

main(int argc, char* argv[])
{
	double *a, *b, *local_a, *local_b;
	double product = 0.0, 
	double local_product = 0.0;

	double start, stop;

	int i;
	int proc_rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0) {
		printf("Enter the number of elements in each vector: "); scanf("%d", &n);
		a = (double *)calloc(N, sizeof(double));
		b = (double *)calloc(N, sizeof(double));
		
		for (i = 0; i < N; i++)
		{
			a[i] = rand() % 10;
			b[i] = rand() % 10;
		}

		start = MPI_Wtime();
	}

	const int block_size = N / size;

	local_a = (double *)calloc(n_bar, sizeof(double));
	local_b = (double *)calloc(n_bar, sizeof(double));

	MPI_Scatter(a, block_size, MPI_DOUBLE, local_a, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, block_size, MPI_DOUBLE, local_b, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i < block_size; i++)
		local_product += local_a[i] * local_b[i];
	free(local_a); 
	free(local_b);

	MPI_Reduce(&local_product, &product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (proc_rank == 0) {
		stop = MPI_Wtime();
		printf("Parallel dot product of %d elements: %f\n", N, product);
		printf("Time: %lf\n", stop - start);

	}

	if (proc_rank == 0) {
		start = MPI_Wtime();
		double seq_product = 0.0;
		for (i = 0; i < N; i++)
			seq_product += a[i] * b[i];

		stop = MPI_Wtime();
		printf("Seq dot product of %d elements: %f\n", N, seq_product);
		printf("Time: %f\n\n", stop - start);

		free(a);
		free(b);
	}

	MPI_Finalize();
}
