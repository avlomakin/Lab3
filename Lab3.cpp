#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include<iomanip>
#include <random>
#include <mpi.h>
using namespace std;


void random_fill(int* matrix, int size)
{
	for (int i = 0; i < size; i++)
	{
		matrix[i] = ((rand() % 50));
	}
}

void multiply_two_matrices(int* a, int* b, int* c, int l, int m, int n)
{
	for (int i = 0; i < l; i++)
	{
		for (int j = 0; j < n; j++)
		{
			c[i * n + j] = 0;
			for (int k = 0; k < m; k++)
			{
				c[i * n + j] += a[i * m + k] * b[k * n + j];
			}
		}
	}
}

int main(int argc, char **argv)
{
	srand(10);
	int size = 2000;
	int size_sq = size * size;
		
	int comm_size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int chunk = size / comm_size;
	int elements_per_proc_count = chunk * size;

	int* a = (int*)malloc(sizeof(int) * size_sq);
	int* b = (int *)malloc(sizeof(int) * size_sq);

	double start = 0;
	
	if (rank == 0)
	{
		start = MPI_Wtime();

		random_fill(a, size_sq);
		random_fill(b, size_sq);
	}

	MPI_Bcast(b, size_sq, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Win win;
	MPI_Win_create(a, sizeof(int) * size_sq, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
	MPI_Win_fence(0, win);

	int* a_proc_chunk = (int*)malloc(sizeof(int) * elements_per_proc_count);
	MPI_Get(a_proc_chunk, elements_per_proc_count, MPI_INT, 0, rank * elements_per_proc_count, elements_per_proc_count, MPI_INT, win);
	MPI_Win_fence(0, win);

	int* result_per_proc = (int*)malloc(sizeof(int) * elements_per_proc_count);
	multiply_two_matrices(a_proc_chunk, b, result_per_proc, chunk, size, size);

	int* result = (int*)malloc(sizeof(int) * size_sq);
	MPI_Gather(result_per_proc, elements_per_proc_count, MPI_INT, result, elements_per_proc_count, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0)
	{
		double endTime = MPI_Wtime();
		double time = endTime - start;
		printf("calculation is finished in %.6f\n", time);

		free(a);
		free(b);
		free(result);
	}

	MPI_Finalize();

	free(a_proc_chunk);
	free(result_per_proc);

	return EXIT_SUCCESS;
}
