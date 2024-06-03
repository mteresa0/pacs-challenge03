INC_DIR = inc
SRC_DIR = src
INCLUDE = -I$(INC_DIR)

CXX = mpicxx
CPPFLAGS = -fopenmp -Wall -Wextra $(INCLUDE)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXEC = main

nodes = 11 21 31 51
threads = 2 3 4 6
process = 2 3 4 6

run : main
	mpirun -np 2 env OMP_NUM_THREAD=6 ./main $(nodes)

all : clean
	make run

main: $(SRCS)
	$(CXX) $(CPPFLAGS) $^ -o main

.PHONY : all clean

clean : 
	rm -f main