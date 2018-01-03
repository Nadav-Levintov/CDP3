#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

/* Debug prints */
#define PRINT(msg) 	printf("%s\n", msg); \
								fflush(stdout)

#define PRINT1(msg,A) 	printf("%s: %d\n", msg,A); \
								fflush(stdout)

#define PRINT2(msg,A,B) 	printf("%s: %d, %d\n", msg,A,B); \
								fflush(stdout)

#define PRINT_ARRAY(arr, len) 	for(int index_bla = 0; index_bla < len; index_bla++) { \
													printf("%d", arr[index_bla]); \
													fflush(stdout); \
												} \
												PRINT("")

#define CALC_DISTANCE(x1, y1, x2, y2) (1+(2*abs(x1-x2))+(2*abs(y1-y2)))
#define MPI_EXEC(rc) if (rc != MPI_SUCCESS) { \
                    PRINT("MPI_Error!"); \
                    MPI_Abort(MPI_COMM_WORLD, rc); }

#define MY_TAG 555
#define INIT_STRUCT struct { \
    int prefix_length; \
    int num_of_prefixes; \
    int xCoord[citiesNum]; \
    int yCoord[citiesNum]; \
    int first_prefix[prefix_length]; \
    }
#define INIT_STRUCT_ARR_SIZE 5

#define RES_STRUCT struct{\
                            int cost;\
                            int path[citiesNum];}


#define RES_STRUCT_ARR_SIZE 2


// Declarations
// ------------
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]);
void recv_and_do_work();
int recv_cities_num();
int calc_initial_lower_bound(int* prefix, int prefix_length, int citiesNum, int **adj_matrix);
int scatter_and_gather(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]);
void
scatter_data(int citiesNum, int xCoord[], int yCoord[], int *prefix_root, int *res_length, int *num_of_root_prefix);
void send_cities_num(int citiesNum, int processesNum);
int calc_length_of_prefix(int citiesNum, int num_of_processes);
int calc_prefix_amount(int citiesNum, int prefix_length);
void build_derived_type(MPI_Datatype *my_data, int citiesNum, int prefix_length);
void next_prefix(int prefix[], int prefix_length, int citiesNum, int step_size);
int next_city(int cityNum, int *prefixCell, int currentCity, bool *visited);
void do_work(void *data_ptr, int citiesNum, int prefix_length, int *final_best_cost, int *final_best_path);
int find_min(int citiesNum, int **adj, int curr_index);
int find_second_min(int citiesNum, int **adj, int curr_index);

// Utils
// ------------

int find_min(int citiesNum, int **adj, int curr_index) {
	int min = INT_MAX;
	for (int k = 0; k < citiesNum; k++)
		if (adj[curr_index][k] < min && curr_index != k)
			min = adj[curr_index][k];

	return min;
}

int find_second_min(int citiesNum, int **adj, int curr_index) {
	int first = INT_MAX, second = INT_MAX;
	for (int j = 0; j < citiesNum; j++) {
		if (curr_index == j)
			continue;
		if (adj[curr_index][j] <= first) {
			second = first;
			first = adj[curr_index][j];
		}
		else if (adj[curr_index][j] <= second &&
			adj[curr_index][j] != first)
			second = adj[curr_index][j];
	}

	return second;
}

void create_adj_matrix(int citiesNum, int xCoord[], int yCoord[], int ***adj_matrix) {
	*adj_matrix = malloc(citiesNum * sizeof(int *));
	for (int i = 0; i < citiesNum; i++) {
		(*adj_matrix)[i] = malloc(citiesNum * sizeof(int));
		for (int j = 0; j < citiesNum; j++) {
			(*adj_matrix)[i][j] = (i == j) ? 0 : CALC_DISTANCE(xCoord[i], yCoord[i], xCoord[j], yCoord[j]);
		}
	}
}

/*
Params:
citiesNum - number of cities
current_path - pointer to array containing the indexes of the cities we already visited by order.
current_index - current index in the path array.
bound - current lower bound.
initial_cost - cost of the previous iteration.
adv_matrix - adj_matrix...
final_res - pointer to int, will hold the cost of the best path found.
final_path - pointer to array that will hold the best path found.
*/
void
branch_and_bound(int citiesNum, int *current_path, int current_index, int bound, int initial_cost, int **adj_matrix,
	int *final_res, int *final_path) {

	int lower_bound = bound, current_cost = initial_cost;

	bool *visited = malloc(citiesNum * sizeof(bool));
	memset(visited, false, sizeof(bool) * citiesNum);
	for (int i = 0; i < current_index; i++) {
		visited[current_path[i]] = true;
	}

	if (current_index == citiesNum) {
		current_cost += adj_matrix[current_path[current_index - 1]][current_path[0]];
		if (current_cost < *final_res) {
			*final_res = current_cost;

			memcpy(final_path, current_path, citiesNum * sizeof(int));
		}
		return;
	}

	for (int i = 0; i < citiesNum; i++) {
		if (visited[i] == false) {
			int temp = lower_bound;
			current_cost += adj_matrix[current_path[current_index - 1]][i];
			/* This is not relevent for us because we will call this function with a prefix, we
			need to use the first case when we build the prefix */
			if (current_index == 1) {
				lower_bound -= ((find_min(citiesNum, adj_matrix, current_path[current_index - 1]) +
					find_min(citiesNum, adj_matrix, i)) / 2);
			}
			else {
				lower_bound -= ((find_second_min(citiesNum, adj_matrix, current_path[current_index - 1]) +
					find_min(citiesNum, adj_matrix, i)) / 2);
			}

			/* This is where we prune, if our current final res is better than the current cost we will not continue */
			if (lower_bound + current_cost < *final_res) {
				current_path[current_index] = i;
				visited[i] = true;

				branch_and_bound(citiesNum, current_path, current_index + 1, lower_bound, current_cost,
					adj_matrix, final_res, final_path);
			}

			/* Continue traversing the tree */
			current_cost -= adj_matrix[current_path[current_index - 1]][i];
			lower_bound = temp;

			memset(visited, false, sizeof(bool) * citiesNum);
			for (int j = 0; j < current_index; j++)
				visited[current_path[j]] = true;
		}
	}
}

// The static parallel algorithm main function.
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]) {
	int myRank;
	MPI_EXEC(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
	if (myRank == 0) {
		// distribute the data, do work, and collect the result
		int result = scatter_and_gather(citiesNum, xCoord, yCoord, shortestPath);
		return result;
	}
	else {
		// calculate partial result and send to CPU 0
		recv_and_do_work();
	}
	return -1;
}

int scatter_and_gather(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]) {
	int prefix[citiesNum];
	int prefix_length = 0, num_of_prefixes = 0;
	int result = 0;
	memset(prefix, 0, sizeof(int) * citiesNum);
	scatter_data(citiesNum, xCoord, yCoord, prefix, &prefix_length, &num_of_prefixes);

	typedef INIT_STRUCT Data;
	Data data;
	data.prefix_length = prefix_length;
	data.num_of_prefixes = num_of_prefixes;
	memmove(data.xCoord, xCoord, sizeof(int) * citiesNum);
	memmove(data.yCoord, yCoord, sizeof(int) * citiesNum);
	memmove(data.first_prefix, prefix, sizeof(int) * prefix_length);
	do_work(&data, citiesNum, prefix_length, &result, shortestPath);
	return result;
}

void
scatter_data(int citiesNum, int xCoord[], int yCoord[], int *prefix_root, int *res_length, int *num_of_root_prefix) {
	int num_of_processes = 0, size = 0;
	MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
	send_cities_num(citiesNum, num_of_processes);
	int prefix_length = calc_length_of_prefix(citiesNum, num_of_processes);
	int num_of_prefixes = calc_prefix_amount(citiesNum, prefix_length);
	typedef INIT_STRUCT Data;
	MPI_Datatype my_data;
	build_derived_type(&my_data, citiesNum, prefix_length);

	int prefix[prefix_length];
	for (int i = 0; i < prefix_length; ++i) {
		prefix[i] = i;
	}
	MPI_EXEC(MPI_Pack_size(1, my_data, MPI_COMM_WORLD, &size));
	size += MPI_BSEND_OVERHEAD;
	void *buf = malloc(size);
	for (int i = 1; i < num_of_processes; ++i) {
		Data data;
		data.num_of_prefixes = num_of_prefixes / num_of_processes;
		if (i <= num_of_prefixes % num_of_processes) {
			data.num_of_prefixes++;
		}
		data.prefix_length = prefix_length;
		memmove(data.xCoord, xCoord, sizeof(int) * citiesNum);
		memmove(data.yCoord, yCoord, sizeof(int) * citiesNum);
		memmove(data.first_prefix, prefix, sizeof(int) * prefix_length);
		MPI_EXEC(MPI_Buffer_attach(buf, size));
		MPI_EXEC(MPI_Bsend(&data, 1, my_data, i, MY_TAG, MPI_COMM_WORLD));
		MPI_EXEC(MPI_Buffer_detach(&buf, &size));
		next_prefix(prefix, prefix_length, citiesNum, data.num_of_prefixes);
	}

	memmove(prefix_root, prefix, sizeof(int) * prefix_length);
	*num_of_root_prefix = num_of_prefixes / num_of_processes;
	*res_length = prefix_length;
	free(buf);
}

void next_prefix(int prefix[], int prefix_length, int citiesNum, int step_size) {
	int curr = prefix_length - 1, advanced = 0;
	bool *visited = malloc(sizeof(bool) * (citiesNum + 1));
	for (int i = 0; i < citiesNum; i++) {
		visited[i] = false;
	}
	for (int i = 0; i < prefix_length; i++) {
		visited[prefix[i]] = true;
	}

	while (step_size > advanced) {
		while (next_city(citiesNum, (prefix + curr), prefix[curr], visited) == citiesNum) {
			curr--;
		}
		while (prefix_length - 1 > curr) {
			curr++;
			next_city(citiesNum, (prefix + curr), prefix[curr], visited);
		}
		advanced++;
	}
	free(visited);
}

int next_city(int citiesNum, int *prefix, int curr_city, bool *visited) {
	visited[curr_city] = false;
	curr_city = (curr_city + 1) % (citiesNum + 1);
	while (curr_city != citiesNum) {
		if (!visited[curr_city]) {
			visited[curr_city] = true;
			break;
		}
		curr_city++;
		curr_city %= (citiesNum + 1);
	}
	*prefix = curr_city;
	return curr_city;
}

void build_derived_type(MPI_Datatype *my_data, int citiesNum, int prefix_length) {
	typedef INIT_STRUCT Data;
	Data data;

	int lengths[INIT_STRUCT_ARR_SIZE];
	MPI_Aint disp[INIT_STRUCT_ARR_SIZE];
	MPI_Datatype types[INIT_STRUCT_ARR_SIZE];
	MPI_Aint addrs[INIT_STRUCT_ARR_SIZE + 1];

	for (int i = 0; i < INIT_STRUCT_ARR_SIZE; ++i) {
		types[i] = MPI_INT;
	}

	lengths[0] = lengths[1] = 1;
	lengths[2] = lengths[3] = citiesNum;
	lengths[4] = prefix_length;

	MPI_EXEC(MPI_Get_address(&data, &addrs[0]));
	MPI_EXEC(MPI_Get_address(&(data.prefix_length), &addrs[1]));
	MPI_EXEC(MPI_Get_address(&(data.num_of_prefixes), &addrs[2]));
	MPI_EXEC(MPI_Get_address(&(data.xCoord), &addrs[3]));
	MPI_EXEC(MPI_Get_address(&(data.yCoord), &addrs[4]));
	MPI_EXEC(MPI_Get_address(&(data.first_prefix), &addrs[5]));

	for (int j = 0; j < INIT_STRUCT_ARR_SIZE; ++j) {
		disp[j] = addrs[j + 1] - addrs[0];
	}

	MPI_EXEC(MPI_Type_create_struct(INIT_STRUCT_ARR_SIZE, lengths, disp, types, my_data));
	MPI_EXEC(MPI_Type_commit(my_data));
}

void build_result_type(MPI_Datatype *message_type_ptr, int citiesNum) {
	typedef RES_STRUCT processResult;
	processResult res_data;

	int lengths[RES_STRUCT_ARR_SIZE];
	MPI_Aint disp[RES_STRUCT_ARR_SIZE];
	MPI_Datatype types[RES_STRUCT_ARR_SIZE];
	MPI_Aint addrs[RES_STRUCT_ARR_SIZE + 1];

	types[0] = types[1] = MPI_INT;

	lengths[0] = 1;
	lengths[1] = citiesNum;

	MPI_EXEC(MPI_Get_address(&res_data, &addrs[0]));
	MPI_EXEC(MPI_Get_address(&(res_data.cost), &addrs[1]));
	MPI_EXEC(MPI_Get_address(&(res_data.path), &addrs[2]));

	disp[0] = addrs[1] - addrs[0];
	disp[1] = addrs[2] - addrs[0];

	MPI_EXEC(MPI_Type_create_struct(RES_STRUCT_ARR_SIZE, lengths, disp, types, message_type_ptr));

	MPI_EXEC(MPI_Type_commit(message_type_ptr));
}

int calc_prefix_amount(int citiesNum, int prefix_length) {
	int result = citiesNum - 1;
	prefix_length -= 2;
	citiesNum -= 2;
	while (prefix_length > 0) {
		result *= citiesNum;
		prefix_length--;
		citiesNum--;
	}
	return result;
}

int calc_length_of_prefix(int citiesNum, int num_of_processes) {
	int result = 2, cities_count = citiesNum - 1;
	int prefix_count = cities_count;
	while (num_of_processes > prefix_count) {
		if (result == citiesNum - 1) {
			return result;
		}
		result++;
		cities_count--;
		prefix_count *= cities_count;
	}
	return result;
}

void send_cities_num(int citiesNum, int processesNum) {

	int pack_size = 0;
	MPI_EXEC(MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &pack_size));
	pack_size += MPI_BSEND_OVERHEAD;
	void *buf = malloc(pack_size);
	for (int i = 1; i < processesNum; ++i) {
		MPI_EXEC(MPI_Buffer_attach(buf, pack_size));
		MPI_EXEC(MPI_Bsend(&citiesNum, 1, MPI_INT, i, MY_TAG, MPI_COMM_WORLD));
		MPI_EXEC(MPI_Buffer_detach(&buf, &pack_size));
	}
	free(buf);
}

int recv_cities_num() {
	int citiesNum;
	MPI_Status sts;
	MPI_EXEC(MPI_Recv(&citiesNum, 1, MPI_INT, 0, MY_TAG, MPI_COMM_WORLD, &sts));
	return citiesNum;
}

int calc_initial_cost(int prefix_length, int** adj_mat, int* prefix) {
	int cost = 0;
	for (int i = 0; i < prefix_length - 1; i++)
	{
		cost += adj_mat[prefix[i]][prefix[i + 1]];
	}
	return cost;
}

void do_work(void *data_ptr, int citiesNum, int prefix_length, int *final_best_cost, int *final_best_path) {
	typedef INIT_STRUCT Data;
	Data *data = (Data *)data_ptr;

	int **adj_mat;
	create_adj_matrix(citiesNum, data->xCoord, data->yCoord, &adj_mat);
	MPI_EXEC(MPI_Barrier(MPI_COMM_WORLD));
	int best_cost = INT_MAX;
	int best_path[citiesNum];
	memset(best_path, 0, citiesNum * sizeof(int));

	for (int i = 0; i < data->num_of_prefixes; i++) {
		int res_cost = INT_MAX;
		int current_cost = calc_initial_cost(prefix_length, adj_mat, data->first_prefix);
		int current_bound = calc_initial_lower_bound(data->first_prefix, prefix_length, citiesNum, adj_mat);
		int current_path[citiesNum];
		memset(current_path, 0, sizeof(int) * citiesNum);
		memmove(current_path, data->first_prefix, sizeof(int) * prefix_length);
		int current_result[citiesNum];
		memset(current_result, 0, sizeof(int) * citiesNum);

		branch_and_bound(citiesNum, current_path, data->prefix_length, current_bound, current_cost, adj_mat, &res_cost,
			current_result);

		if (res_cost < best_cost) {
			best_cost = res_cost;

			memmove(best_path, current_result, sizeof(int)*citiesNum);
		}
		next_prefix(data->first_prefix, prefix_length, citiesNum, 1);
	}
	typedef RES_STRUCT Result;
	Result res;
	res.cost = best_cost;
	memmove(res.path, best_path, sizeof(int) * citiesNum);
	MPI_Datatype my_res;
	build_result_type(&my_res, citiesNum);
	int rank = 0, num_of_processes = 0;
	MPI_EXEC(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
	MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
	Result *resultsArray = NULL;
	if (rank == 0) {
		resultsArray = (Result *)malloc(sizeof(Result) * (num_of_processes));
	}

	MPI_EXEC(MPI_Gather(&res, 1, my_res, resultsArray, 1, my_res, 0, MPI_COMM_WORLD));

	if (rank == 0) {
		int bestIndex = 0;
		for (int i = 0; i < num_of_processes; i++) {
			if (resultsArray[i].cost < resultsArray[bestIndex].cost) {
				bestIndex = i;
			}
		}

		*final_best_cost = resultsArray[bestIndex].cost;
		memmove(final_best_path, resultsArray[bestIndex].path, citiesNum * sizeof(int));
		free(resultsArray);
	}
}


void recv_and_do_work() {
	int num_of_processes = 0;
	MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
	int citiesNum = recv_cities_num();
	int prefix_length = calc_length_of_prefix(citiesNum, num_of_processes);
	typedef INIT_STRUCT Data;
	Data data;
	MPI_Datatype my_type;
	build_derived_type(&my_type, citiesNum, prefix_length);
	MPI_Status sts;
	MPI_EXEC(MPI_Recv(&data, 1, my_type, 0, MY_TAG, MPI_COMM_WORLD, &sts));
	do_work(&data, citiesNum, prefix_length, NULL, NULL);
}

int calc_initial_lower_bound(int* prefix, int prefix_length, int citiesNum, int **adj_matrix) {
	int bound = 0;
	for (int i = 0; i < citiesNum; i++) {
		bound += (find_min(citiesNum, adj_matrix, i) + find_second_min(citiesNum, adj_matrix, i));
	}
	bound = (bound & 1) ? (bound / 2) + 1 : bound / 2;
	for (int i = 1; i < prefix_length; i++)
	{
		int amount;
		if (i == 1) {
			amount = ((find_min(citiesNum, adj_matrix, prefix[i - 1]) +
				find_min(citiesNum, adj_matrix, prefix[i])) / 2);
			bound -= amount;
		}
		else {
			amount = ((find_second_min(citiesNum, adj_matrix, prefix[i - 1]) +
				find_min(citiesNum, adj_matrix, prefix[i])) / 2);
			bound -= amount;
		}
	}
	return bound;
}
