INC_DIR = inc
SRC_DIR = src
INCLUDE = -I$(INC_DIR)
TIMES_FILENAME = outputs/times_
PRINT_FILENAME = outputs/printed_data.txt

CXX = mpicxx
CPPFLAGS = -fopenmp -O3 -Wall -Wextra $(INCLUDE)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXEC = main

MAX ?= 13

nodes = 16 32 64 128
threads = 1 2 3 4 6
process = 1 2 3 4 5 6

1p_multithread : init_times_files main
	@for thread in $(threads); do\
		mpirun -np 1 env OMP_NUM_THREADS=$$thread ./main $(TIMES_FILENAME) $(nodes) >> $(TIMES_FILENAME); done

1t_multiprocess : init_times_files main
	@for proc in $(process); do\
		mpirun -np $$proc env OMP_NUM_THREADS=1 ./main $(TIMES_FILENAME) $(nodes) >> $(PRINT_FILENAME); done


2t_multiprocess : init_times_files main
	@for proc in $(process); do\
		mpirun -np $$proc env OMP_NUM_THREADS=2 ./main $(TIMES_FILENAME) $(nodes) >> $(PRINT_FILENAME); done

all: clean init_times_files main
	@for proc in $(process); do \
		for thr in $(threads); do \
			product=$$((proc*thr)); \
			if [ $$product -lt $(MAX) ]; then\
				mpirun -np $$proc env OMP_NUM_THREADS=$$thr ./main $(TIMES_FILENAME) $(nodes) >> $(PRINT_FILENAME); \
			fi;\
		done;\
	done

main: $(SRCS)
	$(CXX) $(CPPFLAGS) $^ -o main

clean : clean_outputs
	rm -f main

clean_outputs : 
	rm -f outputs/*.csv
	rm -f outputs/*.vtk
	rm -f $(TIMES_FILENAME)*
	rm -f $(PRINT_FILENAME)

init_times_files : clean_outputs
	@for n in $(nodes); do \
		echo n_procs, n_threads, time\ > $(TIMES_FILENAME)$$n.csv ;\
	done

