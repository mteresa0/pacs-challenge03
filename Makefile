INC_DIR = inc
SRC_DIR = src
INCLUDE = -I$(INC_DIR)
TIMES_FILENAME = outputs/times

CXX = mpicxx
CPPFLAGS = -fopenmp -O3 -Wall -Wextra $(INCLUDE)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXEC = main

nodes = 16 32 64 128
threads = 2 3 4 6
process = 6 4 3 2 1

1p_multithread : main
	@for thread in $(threads); do\
		echo ONE PROCESS - MULTITHREADS > $(TIMES_FILENAME)_1p_multithread_$$thread.txt;\
		mpirun -np 1 env OMP_NUM_THREADS=$$thread ./main $(nodes) >> $(TIMES_FILENAME)_1p_multithread_$$thread.txt; done

1t_multiprocess : main
	@for proc in $(process); do\
		echo ONE PROCESS - MULTITHREADS > $(TIMES_FILENAME)_1t_multiproc_$$proc.txt;\
		mpirun -np $$proc env OMP_NUM_THREADS=1 ./main $(nodes) >> $(TIMES_FILENAME)_1t_multiproc_$$proc.txt; done


2t_multiprocess : main
	@for proc in $(process); do\
		echo ONE PROCESS - MULTITHREADS > $(TIMES_FILENAME)_2t_multiproc_$$proc.txt;\
		mpirun -np $$proc env OMP_NUM_THREADS=2 ./main $(nodes) >> $(TIMES_FILENAME)_2t_multiproc_$$proc.txt; done

total_filename = $(TIMES_FILENAME)_total.txt
all: clean main
	@for proc in $(process); do \
		for thr in $(threads); do \
			mpirun -np $$proc env OMP_NUM_THREADS=$$thr ./main $(nodes) >> $(total_filename); \
		done;\
	done

main: $(SRCS)
	$(CXX) $(CPPFLAGS) $^ -o main

clean : 
	rm -f main
	rm -f outputs/*.csv
	rm -f outputs/*.vtk
	rm -f $(TIMES_FILENAME)*