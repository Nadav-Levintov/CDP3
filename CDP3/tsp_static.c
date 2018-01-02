#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

#define CALC_DISTANCE(x1, y1, x2, y2) (1+(2*abs(x1-x2))+(2*abs(y1-y2)))
#define MPI_EXEC(rc) if (rc != MPI_SUCCESS) { \
                    printf("MPI_Error!\n"); \
                    MPI_Abort(MPI_COMM_WORLD, rc); }
#define MY_TAG 555
#define INIT_STRUCT struct { \
    int citiesNum; \
    int prefix_length; \
    int num_of_prefixes; \
    int xCoord[citiesNum]; \
    int yCoord[citiesNum]; \
    int first_prefix[prefix_length]; \
};

// Declarations
// ------------
void receiveAndDoWork();
int calc_initial_lower_bound(int citiesNum, int xCoord[], int yCoord[]);
int distributeAndCollect(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]);
void
distributeData(int citiesNum, int xCoord[], int yCoord[], int *prefix_root, int *res_length, int *num_of_root_prefix);

void sendCitiesNum(int citiesNum, int processesNum);

int calcPrefixLength(int citiesNum, int num_of_processes);

int calcPrefixNum(int citiesNum, int prefix_length);

void build_derived_type(MPI_Datatype *my_data, int citiesNum, int prefix_length);

void next_prefix(int prefix[], int prefix_length, int citiesNum, int step_size);

// Utils
// ------------
int firstMin(int citiesNum, int **adj, int i) {
    int min = INT_MAX;
    for (int k = 0; k < citiesNum; k++)
        if (adj[i][k] < min && i != k)
            min = adj[i][k];
    return min;
}

int secondMin(int citiesNum, int **adj, int i) {
    int first = INT_MAX, second = INT_MAX;
    for (int j = 0; j < citiesNum; j++) {
        if (i == j)
            continue;
        if (adj[i][j] <= first) {
            second = first;
            first = adj[i][j];
        } else if (adj[i][j] <= second &&
                   adj[i][j] != first)
            second = adj[i][j];
    }
    return second;
}

void create_adj_matrix(int citiesNum, int xCoord[], int yCoord[], int ***adj_matrix) {
    *adj_matrix = malloc(citiesNum * sizeof(int *));
    for (int i = 0; i < citiesNum; i++) {
        *adj_matrix[i] = malloc(citiesNum * sizeof(int));
        for (int j = 0; j < citiesNum; j++) {
            *adj_matrix[i][j] = (i == j) ? 0 : CALC_DISTANCE(xCoord[i], yCoord[i], xCoord[j], yCoord[j]);
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
        assert(current_path[current_index - 1] != current_path[0]);
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
            /* This if is not relevent for us because we will call this function with a prefix, we
                need to use the first case when we build the prefix */
            if (current_index == 1) {
                lower_bound -= ((firstMin(citiesNum, adj_matrix, current_path[current_index - 1]) +
                                 firstMin(citiesNum, adj_matrix, i)) / 2);
            } else {
                lower_bound -= ((secondMin(citiesNum, adj_matrix, current_path[current_index - 1]) +
                                 firstMin(citiesNum, adj_matrix, i)) / 2);
            }

            /* This is where we prune, if out current final res is better than the current cost we will not continue */
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
    int rc, myRank;
    MPI_EXEC(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    if (myRank == 0) {
        // distribute the data, do work, and collect the result
        int result = distributeAndCollect(citiesNum, xCoord, yCoord, shortestPath);
        return result;
    } else {
        // calculate partial result and send to CPU 0
        receiveAndDoWork();
    }
    return -1;
}

int distributeAndCollect(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]) {
    int prefix[citiesNum];
    memset(prefix, 0, sizeof(int) * citiesNum);
    // int** adj_mat;
    // create_adj_matrix(citiesNum, xCoord, yCoord, &adj_mat);
    distributeData(citiesNum, xCoord, yCoord, NULL, NULL, NULL);
    Data data;
    calcFirstData(&data, citiesNum); // maybe inside distributeData?
    doWork(data);
    int result = 0;
    CollectResult(shortestPath, &result);
    return result;
}

void
distributeData(int citiesNum, int xCoord[], int yCoord[], int *prefix_root, int *res_length, int *num_of_root_prefix) {
    int num_of_processes = 0, size = 0;
    MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
    sendCitiesNum(citiesNum, num_of_processes);
    int prefix_length = calcPrefixLength(citiesNum, num_of_processes);
    int num_of_prefixes = calcPrefixNum(citiesNum, prefix_length);
    typedef INIT_STRUCT Data;
    MPI_Datatype my_data;
    build_derived_type(&my_data, citiesNum, prefix_length);

    int prefix[prefix_length];
    for (int i = 0; i < prefix_length; ++i) {
        prefix[i] = i;
    }
    MPI_EXEC(MPI_Pack_size(1, my_data, MPI_COMM_WORLD, &size));
    size += MPI_BSEND_OVERHEAD;
    void* buf = malloc(size);
    for (int i = 1; i < num_of_processes; ++i) {
        Data data;
        data.citiesNum = citiesNum;
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
    memmove(prefix_root, prefix, prefix_length);
    *num_of_root_prefix = num_of_prefixes / num_of_processes;
    *res_length = prefix_length;
    free(buf);
}

// TODO: refactor this shit
void next_prefix(int prefix[], int prefix_length, int citiesNum, int step_size) {
    int currentCell = prefixLength - 1;
    int alreadyAdvanced = 0;
    int *cityMapping=(int*)malloc(sizeof(int)*(city_amount + 1));
    for(int i=0; i<city_amount; i++){
        cityMapping[i]=0;
    }
    for(int i=0; i<prefixLength; i++){
        cityMapping[prefix[i]]=1;
    }
    while(alreadyAdvanced < advance_amount){
        while(getNextCity(cityMapping,prefix[currentCell],city_amount,(prefix+currentCell)) == city_amount)
        {
            currentCell--;

        }
        while(currentCell<prefixLength-1)
        {
            currentCell++;
            getNextCity(cityMapping,prefix[currentCell],city_amount,(prefix+currentCell));
        }
        alreadyAdvanced++;
    }
    free(cityMapping);
}

void build_derived_type(MPI_Datatype *my_data, int citiesNum, int prefix_length) {
    typedef INIT_STRUCT Data;
    Data data;

    int block_length[6];
    MPI_Aint displacements[6];
    MPI_Datatype typelist[6];
    MPI_Aint addresses[7];

    for (int i = 0; i < 6; ++i) {
        typelist[i] = MPI_INT;
    }

    block_length[0] = block_length[1] = block_length[2] = 1;
    block_length[3] = block_length[4] = citiesNum;
    block_length[5] = prefix_length;

    MPI_EXEC(MPI_Address(&data, &addresses[0]));
    MPI_EXEC(MPI_Address(&(data.citiesNum), &addresses[1]));
    MPI_EXEC(MPI_Address(&(data.prefixSize), &addresses[2]));
    MPI_EXEC(MPI_Address(&(data.prefixAmount), &addresses[3]));
    MPI_EXEC(MPI_Address(&(data.xCoord), &addresses[4]));
    MPI_EXEC(MPI_Address(&(data.yCoord), &addresses[5]));
    MPI_EXEC(MPI_Address(&(data.firstPrefix), &addresses[6]));

    for (int j = 0; j < 6; ++j) {
        displacements[j] = addresses[j+1] - addresses[0];
    }

    MPI_EXEC(MPI_Type_struct(6, block_length, displacements, typelist, my_data));
    MPI_EXEC(MPI_Type_commit(my_data));
}

int calcPrefixNum(int citiesNum, int prefix_length) {
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

int calcPrefixLength(int citiesNum, int num_of_processes) {
    int result = 2, cities_count = citiesNum - 1;
    int prefix_count = cities_count;
    while (num_of_processes > prefix_count) {
        if(result == citiesNum - 1) {
            return result;
        }
        result++;
        cities_count--;
        prefix_count *= cities_count;
    }
    return result;
}

void sendCitiesNum(int citiesNum, int processesNum) {
    int pack_size = 0;
    MPI_EXEC(MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &size));
    size += MPI_BSEND_OVERHEAD;
    void* buf = malloc(size);
    for (int i = 1; i < processesNum; ++i) {
        MPI_EXEC(MPI_Buffer_attach(buf, size));
        MPI_EXEC(MPI_Bsend(&citiesNum, 1, MPI_INT, i, MY_TAG, MPI_COMM_WORLD));
        MPI_EXEC(MPI_Buffer_detach(&buf, &size));
    }
    free(buf);
}

void receiveAndDoWork() {
    int num_of_processes = 0;
    MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
    int citiesNum = getCitiesNum();
    int currIndex = calcPrefixLength(citiesNum, num_of_processes);
    Data data;
    MPI_Datatype my_type;
    build_derived_type(&data, &my_type, citiesNum, currIndex);
    MPI_Status res;
    MPI_EXEC(MPI_Recv(&data, 1, my_type, 0, MY_TAG, MPI_COMM_WORLD, &res));
    doWork(data);
}
