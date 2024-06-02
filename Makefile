INC_DIR = inc
SRC_DIR = src
INCLUDE = -I$(INC_DIR)

CXX = mpicxx
CPPFLAGS = -fopenmp -Wall -Wextra $(INCLUDE)

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXEC = main

run : main
	mpirun -np 4 ./main 20

all : clean
	make run

main: $(SRCS)
	$(CXX) $(CPPFLAGS) $^ -o main

.PHONY : all clean

clean : 
	rm -f main