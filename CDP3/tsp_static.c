#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>

#define CALC_DISTANCE(x1,y1,x2,y2) (1+(2*abs(x1-x2))+(2*abs(y1-y2)))


int firstMin(int citiesNum, int** adj, int i)
{
	int min = INT_MAX;
	for (int k = 0; k < citiesNum; k++)
		if (adj[i][k] < min && i != k)
			min = adj[i][k];
	return min;
}


int secondMin(int citiesNum, int** adj, int i)
{
	int first = INT_MAX, second = INT_MAX;
	for (int j = 0; j < citiesNum; j++)
	{
		if (i == j)
			continue;

		if (adj[i][j] <= first)
		{
			second = first;
			first = adj[i][j];
		}
		else if (adj[i][j] <= second &&
			adj[i][j] != first)
			second = adj[i][j];
	}
	return second;
}

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


/*
	Params:
		citiesNum - number of cities
		current_path - pointer to array containing the indexes of the cities we already visited by order.
		current_index - current index in the path array.
		bound - current lower bound.
		initial_cost - cost of the previeus iteration.
		adv_matrix - adj_matrix...
		final_res - pointer to int, will hold the cost of the best path found.
		final_path - pointer to array that will hold the best path found.


*/
void branch_and_bound(int citiesNum, int *current_path, int current_index, int bound, int initial_cost, int** adj_matrix, int *final_res, int* final_path)
{
	int lower_bound = bound, current_cost = initial_cost;

	bool* visited = malloc(citiesNum * sizeof(bool));
	memset(visited, false, sizeof(bool)*citiesNum);
	for (int i = 0; i < current_index; i++)
	{
		visited[current_path[i]] = true;
	}

	if (current_index == citiesNum)
	{
		assert(current_path[current_index - 1] != current_path[0]);
		current_cost += adj_matrix[current_path[current_index - 1]][current_path[0]];
		if (current_cost < *final_res)
		{
			*final_res = current_cost;
			memcpy(final_path, current_path, citiesNum * sizeof(int));
		}
		return;
	}

	for (int i = 0; i < citiesNum; i++)
	{
		if (visited[i] == false)
		{
			int temp = lower_bound;
			current_cost += adj_matrix[current_path[current_index - 1]][i];
			/* This if is not relevent for us because we will call this function with a prefix, we 
				need to use the first case when we build the prefix */
			if (current_index == 1)
			{
				lower_bound -= ((firstMin(citiesNum, adj_matrix, current_path[current_index - 1]) +
					firstMin(citiesNum, adj_matrix, i)) / 2);
			}
			else
			{
				lower_bound -= ((secondMin(citiesNum, adj_matrix, current_path[current_index - 1]) +
					firstMin(citiesNum, adj_matrix, i)) / 2);
			}

			/* This is where we prune, if out current final res is better than the current cost we will not continue */
			if (lower_bound + current_cost < *final_res)
			{
				current_path[current_index] = i;
				visited[i] = true;

				branch_and_bound(citiesNum, current_path, current_index + 1, lower_bound, current_cost,
					adj_matrix, final_res, final_path);
			}

			/* Continue traversing the tree */
			current_cost -= adj_matrix[current_path[current_index - 1]][i];
			lower_bound = temp;

			memset(visited, false, sizeof(bool)*citiesNum);
			for (int j = 0; j < current_index; j++)
				visited[current_path[j]] = true;

		}
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
		int** adj_matrix;
		create_adj_matrix(citiesNum, xCoord, yCoord, &adj_matrix);
		//divide work
	}
	else {
		//Do work
	}



	return -1;	//TODO
}
