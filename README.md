# ParallelProgrammingExample

This Project contains two examples for Parallel Programming. 
The first example try to simulate the `N-Body` problem by Numerical Method.
The second example try to simulate the `TSP` problem.

Every example includes three implements: `OpenMP`, `MPI` and `CUDA`.

The main program uses `GMP` for high precision numerical calculation.

## N-Body

When you are going to compile this example, you should have installed
`CUDA`, `openmpi`, `openmp` and other tools.

There are some tutorials: [MPI Tutorial](http://1.117.83.246/index.php/archives/266/),
[OpenMP Tutorial](http://1.117.83.246/index.php/archives/265/).

`$ cd n-body`

To compile the code:

`$ cmake -B build`

`$ cd build && make`

To enter the path (`build/src`):

`$ cd src`

### OpenMP

`$ ./nbody_openmp < ../../3_Body_test.txt`

### MPI

`$ ./nbody_mpi < ../../3_Body_test.txt`

`$ mpiexec -n 3 --use-hwthread-cpus ./nbody_mpi < ../../3_Body_test.txt`

### CUDA


### 3_Body_test.txt

This is a testcase contains the mass, initial-position and 
initial-velocity of sun, earth and moon.

## TSP

