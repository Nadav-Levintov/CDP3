#include <mpi.h>
#include <mpi.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>

/* Debug prints */
#define PRINT(msg)    printf("%s\n", msg); \
                                fflush(stdout)

#define PRINT1(msg, A)    printf("%s: %d\n", msg,A); \
                                fflush(stdout)

#define PRINTP(msg, A)    printf("%s: %p\n", msg,A); \
                                fflush(stdout)

#define PRINT2(msg, A, B)    printf("%s: %s, %d\n", msg,A,B); \
                                fflush(stdout)

#define PRINT_ARRAY(arr, len)    for(int index_bla = 0; index_bla < len; index_bla++) { \
                                                    printf("%d", arr[index_bla]); \
                                                    fflush(stdout); \
                                                } \
                                                PRINT("")

/* Other Macros */
#define CALC_DISTANCE(x1, y1, x2, y2) (1+(2*abs(x1-x2))+(2*abs(y1-y2)))
#define MPI_EXEC(rc) if (rc != MPI_SUCCESS) { \
                    PRINT2("MPI_Error! ", __func__, __LINE__); \
                    MPI_Abort(MPI_COMM_WORLD, rc); }

#define MY_TAG 555
#define TAG_KILL 556
#define TAG_WORKER_RES 557
#define TAG_NEW_BEST_COST 558
#define TAG_NEW_JOB 559
#define DYANAMIC_WORKER_FACTOR 100

#define JOB_STRUCT struct{\
                            int is_valid; \
                            int prefix_length;\
                            int prefix[citiesNum];}
#define JOB_STRUCT_ARR_SIZE 3

#define WORKER_RES_STRUCT struct{\
                                int rank;\
                                int res_cost;\
                                int res_path[citiesNum];}
#define WORKER_RES_STRUCT_ARR_SIZE 3

#define COST_STRUCT struct{\
                                int cost;}
#define COST_STRUCT_ARR_SIZE 1

/* function declaretions */

int next_city(int cityNum, int *prefixCell, int currentCity, bool *visited);

void next_prefix(int prefix[], int prefix_length, int citiesNum, int step_size);

int find_min(int citiesNum, int **adj, int curr_index);

int find_second_min(int citiesNum, int **adj, int curr_index);

int calc_initial_lower_bound(int *prefix, int prefix_length, int citiesNum, int **adj_matrix);

void
branch_and_bound(int citiesNum, int *current_path, int current_index, int bound, int initial_cost, int **adj_matrix,
                 int *final_res, int *final_path);

void create_adj_matrix(int citiesNum, int xCoord[], int yCoord[], int ***adj_matrix);

int master_func(int citiesNum, int *shortestPath, int num_of_workers);

void worker_func(int citiesNum, int **adj_matrix, int num_of_processes);

int calc_length_of_prefix(int citiesNum, int num_of_processes);

int calc_prefix_amount(int citiesNum, int prefix_length);

void build_cost_type(MPI_Datatype *cost_type);

int calc_initial_cost(int prefix_length, int **adj_mat, int *prefix);
void build_job_type(MPI_Datatype *message_type_ptr, int citiesNum);
void build_worker_res_type(MPI_Datatype *worker_res_type, int citiesNum);
void free_adj_matrix(int citiesNum, int ***adj_matrix);
// Utils
// ------------



void checkForBetterResult(int *my_lowest_cost);

void sendMyBetterResult(int rank, int num_of_processes, int my_lowest_cost, void *send_buffer, int send_buffer_size,
                        MPI_Datatype cost_type);

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
        } else if (adj[curr_index][j] <= second &&
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

void free_adj_matrix(int citiesNum, int ***adj_matrix) {
    for (int i = 0; i < citiesNum; i++) {
        free((*adj_matrix)[i]);
    }
    free(*adj_matrix);
}

/*
* Branch and Bound function
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
            } else {
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
    free(visited);
}

/* calculate the lower bound given a prefix */
int calc_initial_lower_bound(int *prefix, int prefix_length, int citiesNum, int **adj_matrix) {
    int bound = 0;
    for (int i = 0; i < citiesNum; i++) {
        bound += (find_min(citiesNum, adj_matrix, i) + find_second_min(citiesNum, adj_matrix, i));
    }
    bound = (bound & 1) ? (bound / 2) + 1 : bound / 2;
    for (int i = 1; i < prefix_length; i++) {
        int amount;
        if (i == 1) {
            amount = ((find_min(citiesNum, adj_matrix, prefix[i - 1]) +
                       find_min(citiesNum, adj_matrix, prefix[i])) / 2);
            bound -= amount;
        } else {
            amount = ((find_second_min(citiesNum, adj_matrix, prefix[i - 1]) +
                       find_min(citiesNum, adj_matrix, prefix[i])) / 2);
            bound -= amount;
        }
    }
    return bound;
}


// The dynamic parellel algorithm main function.
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]) {
    int myRank, num_of_processes = 0;
    MPI_EXEC(MPI_Comm_size(MPI_COMM_WORLD, &num_of_processes));
    if (num_of_processes < 2) {
        return -1;
    }

    int **adj_matrix;
    create_adj_matrix(citiesNum, xCoord, yCoord, &adj_matrix);
    MPI_EXEC(MPI_Comm_rank(MPI_COMM_WORLD, &myRank));
    if (myRank == 0) {
        int result = master_func(citiesNum, shortestPath, num_of_processes);
        free_adj_matrix(citiesNum, &adj_matrix);
        return result;
    } else {
        worker_func(citiesNum, adj_matrix, num_of_processes);
    }
    return -1;
}


/* function implementation */

/* gets a prefix and calculates the next prefix according to the given step size */
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

/* gets a city and returns the next unvisited city */
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

int master_func(int citiesNum, int *shortestPath, int num_of_proccesses) {
    int num_of_workers = 0;
    int initial_prefix_length = calc_length_of_prefix(citiesNum, num_of_workers * 1);
    //initial_prefix_length += 4;
    //initial_prefix_length = initial_prefix_length > citiesNum ? citiesNum : initial_prefix_length;
    int num_of_prefixes = calc_prefix_amount(citiesNum, initial_prefix_length);

    int current_prefix[citiesNum];
    memset(current_prefix, 0, sizeof(int) * citiesNum);
    for (int i = 0; i < initial_prefix_length; i++) {
        current_prefix[i] = i;
    }


    typedef JOB_STRUCT Job;
    MPI_Datatype job_type;
    build_job_type(&job_type, citiesNum);
    Job worker_jobs_arr[num_of_proccesses];

    //starting from 1 because 0 is master
    for (int i = 1; i < num_of_proccesses; i++) {
        if (num_of_prefixes == 0) {
            worker_jobs_arr[i].is_valid = 0;
            continue;
        }
        worker_jobs_arr[i].is_valid = 1;
        worker_jobs_arr[i].prefix_length = initial_prefix_length;
        memmove(worker_jobs_arr[i].prefix, current_prefix, sizeof(int) * citiesNum);
        next_prefix(current_prefix, initial_prefix_length, citiesNum, 1);
        num_of_prefixes--;
        num_of_workers++;
    }

    Job dummy_job;
    MPI_EXEC(MPI_Scatter(worker_jobs_arr, 1, job_type, &dummy_job, 1, job_type, 0, MPI_COMM_WORLD));


    int res_cost = INT_MAX;
    int res_path[citiesNum];
    memset(res_path, 0, citiesNum * sizeof(int));
    MPI_Status status;
    typedef WORKER_RES_STRUCT Res;
    Res worker_res;
    MPI_Datatype Res_type;
    build_worker_res_type(&Res_type, citiesNum);
    Job next_job;
    while (num_of_prefixes) {
        MPI_EXEC(MPI_Recv(&worker_res, 1, Res_type, MPI_ANY_SOURCE, TAG_WORKER_RES, MPI_COMM_WORLD, &status));
        if (worker_res.res_cost < res_cost) {
            res_cost = worker_res.res_cost;
            memmove(res_path, worker_res.res_path, citiesNum * sizeof(int));
        }
        next_job.is_valid = 1;
        next_job.prefix_length = initial_prefix_length;
        memmove(next_job.prefix, current_prefix, sizeof(int) * initial_prefix_length);
        MPI_Ssend(&next_job, 1, job_type, worker_res.rank, TAG_NEW_JOB, MPI_COMM_WORLD);
        next_prefix(current_prefix, initial_prefix_length, citiesNum, 1);
        num_of_prefixes--;
    }

    while (num_of_workers > 0) {
        MPI_EXEC(MPI_Recv(&worker_res, 1, Res_type, MPI_ANY_SOURCE, TAG_WORKER_RES, MPI_COMM_WORLD, &status));
        if (worker_res.res_cost < res_cost) {
            res_cost = worker_res.res_cost;
            memmove(res_path, worker_res.res_path, citiesNum * sizeof(int));
        }
        num_of_workers--;
    }

    for (int i = 1; i < num_of_proccesses; i++) {
        MPI_EXEC(MPI_Ssend(&next_job, 1, job_type, i, TAG_KILL, MPI_COMM_WORLD));
    }
    memmove(shortestPath, res_path, sizeof(int) * citiesNum);

    return res_cost;
}

void build_job_type(MPI_Datatype *message_type_ptr, int citiesNum) {
    typedef JOB_STRUCT Job;
    Job innerData;

    int lengths[JOB_STRUCT_ARR_SIZE];
    MPI_Aint disp[JOB_STRUCT_ARR_SIZE];
    MPI_Datatype types[JOB_STRUCT_ARR_SIZE];
    MPI_Aint addrs[JOB_STRUCT_ARR_SIZE + 1];

    types[0] = types[1] = types[2] = MPI_INT;

    lengths[0] = 1;
    lengths[1] = 1;
    lengths[2] = citiesNum;


    MPI_EXEC(MPI_Get_address(&innerData, &addrs[0]));
    MPI_EXEC(MPI_Get_address(&(innerData.is_valid), &addrs[1]));
    MPI_EXEC(MPI_Get_address(&(innerData.prefix_length), &addrs[2]));
    MPI_EXEC(MPI_Get_address(&(innerData.prefix), &addrs[3]));

    for (int i = 0; i < JOB_STRUCT_ARR_SIZE; ++i) {
        disp[i] = addrs[i + 1] - addrs[0];
    }
    MPI_EXEC(MPI_Type_create_struct(JOB_STRUCT_ARR_SIZE, lengths, disp, types, message_type_ptr));

    MPI_EXEC(MPI_Type_commit(message_type_ptr));
}

void build_worker_res_type(MPI_Datatype *worker_res_type, int citiesNum) {
    typedef WORKER_RES_STRUCT Worker_Res;
    Worker_Res res;
    int lengths[WORKER_RES_STRUCT_ARR_SIZE];
    MPI_Aint disp[WORKER_RES_STRUCT_ARR_SIZE];
    MPI_Datatype types[WORKER_RES_STRUCT_ARR_SIZE];
    MPI_Aint addrs[WORKER_RES_STRUCT_ARR_SIZE + 1];

    for (int i = 0; i < JOB_STRUCT_ARR_SIZE; i++) {
        types[i] = MPI_INT;
    }

    lengths[0] = lengths[1] = 1;
    lengths[2] = citiesNum;

    MPI_EXEC(MPI_Get_address(&res, &addrs[0]));
    MPI_EXEC(MPI_Get_address(&(res.rank), &addrs[1]));
    MPI_EXEC(MPI_Get_address(&(res.res_cost), &addrs[2]));
    MPI_EXEC(MPI_Get_address(&(res.res_path), &addrs[3]));

    disp[0] = addrs[1] - addrs[0];
    disp[1] = addrs[2] - addrs[0];
    disp[2] = addrs[3] - addrs[0];

    MPI_EXEC(MPI_Type_create_struct(JOB_STRUCT_ARR_SIZE, lengths, disp, types, worker_res_type));

    MPI_EXEC(MPI_Type_commit(worker_res_type));
}

void worker_func(int citiesNum, int **adj_matrix, int num_of_processes) {
    int rank;
    MPI_EXEC(MPI_Comm_rank(MPI_COMM_WORLD, &rank));

    MPI_Datatype job_type, cost_type, worker_res_type;
    build_job_type(&job_type, citiesNum);
    build_worker_res_type(&worker_res_type, citiesNum);
    build_cost_type(&cost_type);
    typedef JOB_STRUCT Job;
    typedef COST_STRUCT Cost;
    typedef WORKER_RES_STRUCT Res;
    Job job;
    Cost cost;
    Res res;
    int send_buffer_size = 0;
    MPI_EXEC(MPI_Pack_size(1, cost_type, MPI_COMM_WORLD, &send_buffer_size));
    send_buffer_size += MPI_BSEND_OVERHEAD;
    void* send_buffer = malloc(send_buffer_size * num_of_processes);

    // get first job
    MPI_EXEC(MPI_Scatter(NULL, 1, job_type, &job, 1, job_type, 0, MPI_COMM_WORLD));

    int my_lowest_cost = INT_MAX;
    int my_result_cost = INT_MAX;
    int current_bound, current_cost;
    bool process_alive = true;
    int shortest_path[citiesNum];

    while (process_alive) {
        if (job.is_valid) {
            current_bound = calc_initial_lower_bound(job.prefix, job.prefix_length, citiesNum, adj_matrix);
            current_cost = calc_initial_cost(job.prefix_length, adj_matrix, job.prefix);
            int tmp_cost = my_lowest_cost;
            branch_and_bound(citiesNum, job.prefix, job.prefix_length, current_bound, current_cost, adj_matrix,
                             &tmp_cost, shortest_path);

            checkForBetterResult(&my_lowest_cost);

            if(tmp_cost < my_lowest_cost) {
                // found better result than the rest
                my_lowest_cost = tmp_cost;
                sendMyBetterResult(rank, num_of_processes, my_lowest_cost, send_buffer, send_buffer_size, cost_type);
                my_result_cost = my_lowest_cost;
            } else {
                my_result_cost = INT_MAX;
            }

            res.res_cost = my_result_cost;
            if (my_result_cost < INT_MAX) {
                memmove(res.res_path, shortest_path, citiesNum * sizeof(int));
            }
            res.rank = rank;
            MPI_EXEC(MPI_Ssend(&res, 1, worker_res_type, 0, TAG_WORKER_RES, MPI_COMM_WORLD));
        }
        bool received_new_job = false;
        while (!received_new_job) {
            MPI_Status status;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            switch (status.MPI_TAG) {
                case TAG_NEW_BEST_COST:
                    MPI_EXEC(MPI_Recv(&cost, 1, cost_type, MPI_ANY_SOURCE, TAG_NEW_BEST_COST, MPI_COMM_WORLD, &status));
                    if (cost.cost < my_lowest_cost)
                    {
                        my_result_cost = INT_MAX;
                        my_lowest_cost = cost.cost;
                    }
                    break;
                case TAG_NEW_JOB:
                    MPI_EXEC(MPI_Recv(&job, 1, job_type, 0, TAG_NEW_JOB, MPI_COMM_WORLD, &status));
                    received_new_job = true;

                    break;
                case TAG_KILL:
                    MPI_EXEC(MPI_Recv(&job, 1, job_type, 0, TAG_KILL, MPI_COMM_WORLD, &status));
                    received_new_job = true;
                    process_alive = false;
                    break;
                default:
                    MPI_EXEC(MPI_Recv(&cost, 1, cost_type, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status));
                    job.is_valid = 0;
                    break;
            }
        }
    }
    free(send_buffer);
}

void checkForBetterResult(int *my_lowest_cost) {
    MPI_Datatype cost_type;
    build_cost_type(&cost_type);
    typedef COST_STRUCT Cost;
    Cost cost;
    MPI_Status status;
    int is_available_flag;
    MPI_Iprobe(MPI_ANY_SOURCE, TAG_NEW_BEST_COST, MPI_COMM_WORLD, &is_available_flag,&status);
    while (is_available_flag) {
        // somebody sent new best cost
        MPI_EXEC(MPI_Recv(&cost, 1, cost_type, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status));
        if(cost.cost < *my_lowest_cost) {
            *my_lowest_cost = cost.cost;
        }
        is_available_flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_NEW_BEST_COST, MPI_COMM_WORLD, &is_available_flag,&status);
    }
}

void sendMyBetterResult(int rank, int num_of_processes, int my_lowest_cost, void *send_buffer, int send_buffer_size,
                        MPI_Datatype cost_type) {
    typedef COST_STRUCT Cost;
    Cost cost;
    cost.cost = my_lowest_cost;
    Cost* cost_ptr = (Cost*)send_buffer;
    for (int i = 1; i < num_of_processes; ++i) {
        if(i == rank) {
            continue;
        }
        void* curr_buffer = cost_ptr + i;
        MPI_Buffer_detach(&(curr_buffer), &send_buffer_size);
        MPI_EXEC(MPI_Buffer_attach(curr_buffer, send_buffer_size));
        MPI_EXEC(MPI_Bsend(&cost, 1, cost_type, i, TAG_NEW_BEST_COST, MPI_COMM_WORLD));
    }
}

/* calculates the number if prefixes according to the number of cities */
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

/* calculates the length of each prefix according to the number of cities */
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
    if(result > citiesNum) {
        return citiesNum;
    }
    return result;
}

void build_cost_type(MPI_Datatype *cost_type) {
    typedef COST_STRUCT Cost;
    Cost cost;

    int lengths[COST_STRUCT_ARR_SIZE];
    MPI_Aint disp[COST_STRUCT_ARR_SIZE];
    MPI_Datatype types[COST_STRUCT_ARR_SIZE];
    MPI_Aint addrs[COST_STRUCT_ARR_SIZE + 1];
    types[0] = MPI_INT;

    lengths[0] = 1;

    MPI_EXEC(MPI_Get_address(&cost, &addrs[0]));
    MPI_EXEC(MPI_Get_address(&(cost.cost), &addrs[1]));


    disp[0] = addrs[1] - addrs[0];

    MPI_EXEC(MPI_Type_create_struct(COST_STRUCT_ARR_SIZE, lengths, disp, types, cost_type));
    MPI_EXEC(MPI_Type_commit(cost_type));
}

int calc_initial_cost(int prefix_length, int **adj_mat, int *prefix) {
    int cost = 0;
    for (int i = 0; i < prefix_length - 1; i++) {
        cost += adj_mat[prefix[i]][prefix[i + 1]];
    }
    return cost;
}

