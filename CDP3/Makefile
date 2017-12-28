# compile the c files into o's.
tsp.o: tsp.c
	srun mpicc -c -Wall -std=c99 tsp.c

tsp_static.o: tsp_static.c
	srun mpicc -c -Wall -std=c99 tsp_static.c

main.o: main.c
	srun mpicc -c -Wall -std=c99 main.c

# make the executables: course example, dynamic and static.
tsp_236370: main.o tsp_236370.o
	srun mpicc -o tsp_236370 main.o tsp_236370.o

tsp: main.o tsp.o
	srun mpicc -o tsp main.o tsp.o

tsp_static: main.o tsp_static.o
	srun mpicc -o tsp_static main.o tsp_static.o
	
# clean previous builds
clean:
	srun make cleaning

cleaning:
	@rm -f tsp.o main.o tsp_static.o
	@echo "Cleaned"

# run the executables: course example, dynamic and static.
# you may change the number of procs below to use more than 2.
run_236370: tsp_236370
	@echo "Running supplied implementation: tsp_236370" 
	#mpirun -np 2 tsp_236370
	srun -n 8 tsp_236370
	
run_tsp: tsp
	@echo "Running dynamic implementation: tsp"
	srun -n 8 tsp
	
run_tsp_static: tsp_static
	@echo "Running static implementation: tsp_static"
	srun -n 8 tsp_static
