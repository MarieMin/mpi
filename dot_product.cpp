#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <cstdlib>

using namespace std;

#define N 1000

double dot_product(double *a, double *b, int start_idx, int end_idx) {
	int i;
	double product = 0.0;

	for (i = 0; start_idx < end_idx; i++)
	{
		product += a[i] * b[i];
	}
	return product;
}

int main(int argc, char **argv) {

	double *a, *b;
	double product, local_product;

	double start, stop, parallel_dur, seq_dur;
	double max_time, min_time, avg_time;

	int proc_rank, size;
	int n_elements_recieved;

	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	if (proc_rank == 0)
	{
		printf("\nComm size: %d\n", size);


		int index, i;

		a = (double*)malloc(N * sizeof(double));
		b = (double*)malloc(N * sizeof(double));


		for (i = 0; i < N; i++)
		{
			a[i] = rand() % 10;
			b[i] = rand() % 10;
		}

		const int block_size = N / size;

		start = MPI_Wtime();
		if (size > 1) {
			// distributes the portion of array to child processes
			for (i = 1; i < size - 1; i++) {
				index = i * block_size;

				MPI_Send(&block_size,
					1, MPI_INT, i, 0,
					MPI_COMM_WORLD);

				MPI_Send(&a[index],
					block_size,
					MPI_DOUBLE, i, 1,
					MPI_COMM_WORLD);

				MPI_Send(&b[index],
					block_size,
					MPI_DOUBLE, i, 2,
					MPI_COMM_WORLD);
			}

			// last process adds remaining elements 
			index = i * block_size;
			int elements_left = N - index;

			MPI_Send(&elements_left,
				1, MPI_INT,
				i, 0,
				MPI_COMM_WORLD);

			MPI_Send(&a[index],
				elements_left,
				MPI_DOUBLE, i, 1,
				MPI_COMM_WORLD);

			MPI_Send(&b[index],
				elements_left,
				MPI_DOUBLE, i, 2,
				MPI_COMM_WORLD);
		}

		// master portion of scalar product
		product = dot_product(a, b, 0, block_size);

		// collect local products from other processes 
		int tmp;
		for (i = 1; i < size; i++) {
			MPI_Recv(&tmp, 1, MPI_DOUBLE,
				MPI_ANY_SOURCE, 0,
				MPI_COMM_WORLD,
				&status);
			int sender = status.MPI_SOURCE;

			//printf("Curr local prod of proc [%d] is: %d\n", sender, tmp);

			product += tmp;
		}

		stop = MPI_Wtime();

		parallel_dur = stop - start;
		/*compute max, min, and average timing statistics*/
		MPI_Reduce(&parallel_dur, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&parallel_dur, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&parallel_dur, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		printf("Parallel dot product of %d elements: %f\n", N, product);
		avg_time /= size;
		printf("Min: %lf  Max: %lf  Avg:  %lf\n", min_time, max_time, avg_time);
		

		// sequential section
		start = MPI_Wtime();
		double seq_product = dot_product(a, b, 0, N);
		stop = MPI_Wtime();

		printf("Seq dot product of %d elements: %f\n", N, seq_product);
		printf("Time: %f\n\n", stop - start);
		

	}
	else {
		MPI_Recv(&n_elements_recieved,
			1, MPI_DOUBLE, 0, 0,
			MPI_COMM_WORLD,
			&status);

		double *local_a = (double*)malloc(n_elements_recieved * sizeof(double));
		double *local_b = (double*)malloc(n_elements_recieved * sizeof(double));
	
		MPI_Recv(local_a, n_elements_recieved,
			MPI_DOUBLE, 0, 1,
			MPI_COMM_WORLD,
			&status);

		MPI_Recv(local_b, n_elements_recieved,
			MPI_DOUBLE, 0, 2,
			MPI_COMM_WORLD,
			&status);

		local_product = dot_product(local_a, local_b, 0, n_elements_recieved);
		
		// sends the partial min to the root process 
		MPI_Send(&local_product, 1, MPI_DOUBLE,
			0, 0, MPI_COMM_WORLD);
	}


	MPI_Finalize();
	return 0;
}
