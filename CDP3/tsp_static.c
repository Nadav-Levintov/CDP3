#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>

#define CALC_DISTANCE(x1,y1,x2,y2) (1+(2*abs(x1-x2))+(2*abs(y1-y2)))

void create_adj_matrix(int citiesNum, int xCoord[], int yCoord[], int*** adj_matrix)
{
	*adj_matrix = malloc(citiesNum * sizeof(int*));

	for (int i = 0; i < citiesNum; i++)
	{
		*adj_matrix[i] = malloc(citiesNum * sizeof(int));

		for (int j = 0; j < citiesNum; j++)
		{
			*adj_matrix[i][j] = (i == j) ? 0 : CALC_DISTANCE(xCoord[i], yCoord[i], xCoord[j], yCoord[j]);
		}
	}
}

int calc_initial_lower_bound(int citiesNum, int xCoord[], int yCoord[]);

int branch_and_bound(int citiesNum, int xCoord[], int yCoord[], int *current_path, int current_index, int initial_bound)
{
	int lower_bound = initial_bound, current_cost = 0;

	int** adj_matrix;
	create_adj_matrix(citiesNum, xCoord, yCoord, &adj_matrix);

	bool* visited = malloc(citiesNum * sizeof(bool));
	memset(visited, 0, sizeof(visited));
	for (int i = 0; i < current_index; i++)
	{
		visited[current_path[i]] = true;
	}



}

// The static parellel algorithm main function.
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[])
{
	int rc;
	int number_of_threads, myRank;
	rc = MPI_Comm_size(MPI_COMM_WORLD, &number_of_threads);
	assert(rc == MPI_SUCCESS);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	assert(rc == MPI_SUCCESS);

	if (myRank == 0)
	{
		//divide work
	}
	else {
		//Do work
	}



	return -1;	//TODO
}
