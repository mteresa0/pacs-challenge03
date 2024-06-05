INC_DIR = inc
SRC_DIR = src
INCLUDE = -I$(INC_DIR)
TIMES_FILENAME = outputs/times.txt

CXX = mpicxx
CPPFLAGS = -fopenmp -O3 -Wall -Wextra $(INCLUDE)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXEC = main

nodes = 16 32 64 128 256
threads = 2 3 4 6
process = 2 3 4 6

only_threads : main
	echo 1 PROCESS - MULTITHREADS > $(TIMES_FILENAME)
	@for thread in $(threads); do\
		mpirun -np 1 env OMP_NUM_THREAD=$$thread ./main $(nodes)> $(TIMES_FILENAME); done

only_process : main
	@for proc in $(process); do\
		mpirun -np $$proc env OMP_NUM_THREAD=6 ./main $(nodes) > $(TIMES_FILENAME); done

main: $(SRCS)
	$(CXX) $(CPPFLAGS) $^ -o main

clean : 
	rm -f main
	rm outputs/*