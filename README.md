# Thrid challenge - PACS
Thrid challenge of the course Advaced Programming for Scientific Computing, 23/24, Politecnico di Milano.

## Installation 
To clone the repository run:
```
git clone git@github.com:mteresa0/pacs-challenge03.git
```
## Configuration
First, you have to set `MAX` as the maximum (thread*processes), accoing to you personal computer. You could edit the `Makefile`; for example:
```
MAX=12
```

Alternatively, you may set an environmental variable with the same name using the command 
```
export MAX=12
```

## Usage
In the terminal run:
- `make main` to compile;
- `make 1p_multithread` to run with 1 process and 1 to 6 threads;
- `make 1t_multiprocess` to run with 1 thread and 1 to 6 process; 
- `make 2t_multiprocess` to run with 2 thread and 1 to 6 process; 
- `make all` to compile and run the program with 1 to 6 threads and 1 to 6 process ( if thread x processed are less than `MAX`); 
- `make clean_outputs` to all outputs files.

## Files
- `src/main.cpp` - main file
- `inc/boundaries.hpp` - `boundaries` struct declaration
- `inc/domain.hpp` - `domain` struct declaration
- `inc/solver.hpp` - `solver` class declaration
- `inc/writing.hpp` - declarations of methods `write_VTK_2D_domain` and `write_convergence`
- `src/solver.cpp` - `solver` class definitions (also contains the declaration of `check_boundaries_compatibility`, a method of `boundaries`)
- `src/writing.cpp` -  definitions of `write_VTK_2D_domain` and `write_convergence`
- in folder `outputs` there are the outputs of the program and other files, in particular:
    - `hw.info` processors info
    - `times_256.csv` - times for different processes and threads with a grid 256x256
    - `times_N.csv` - times for different processes and threds with a grid N*N (not pushed)
    - `convergence_nProc_nThr.csv` - L2 norms for different grid sizes with nProc number of processes and nThr number of threads (not pushed)
    - `printed_data.txt` - output of each program run (not pushed)



