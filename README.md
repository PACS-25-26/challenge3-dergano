[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/tKSbaXxd)
# challenge3
The third challenge

# Overview
This challenge solves the Laplace equation that models heat diffusion over a square domain, with Dirichlet boundary condition (either homogenous or not). The code uses the Jacobi method, with the improvment of the Block Jacobi extension.

# Dependencies

**muParser**
This library is used to handle the forcing term and the boundary function, passed by the user via command line. The library allows us to evaluate a string, giving us a vector of double. It need to be installed.

**OpenMP and MPI**
Recquired for parallelization.


# Building
The project includes a standard `Makefile` to simplify the build process. To compile the code, open your terminal in the project's root directory and run:

```bash
make
```

This command will compile the source files located in the src/ directory and generate an executable named laplace in the root directory.
To remove the compiled object files (.o) and the laplace executable, run:

```bash
make clean
```

# Running 
Since the program uses MPI for distributed memory parallelism and OpenMP for shared memory, it must be executed using `mpirun` (or `mpiexec`).

```bash
mpirun -np <num_mpi_processes> ./laplace <k_value> <num_threads> "<forcing_function>" "<boundary_function>"
```
# Testing
A bash script is provided to automate scalability testing across different serial, pure MPI, pure OpenMP, and hybrid configurations. To execute it, change the directory into the test folder. Then execute (the first only one time):

```bash
chmod +x run_scalability.sh
./run_scalability.sh
```

The results are collected in the scalability_report file. It's provided also a RESULT file wiith some comments about the results.

# Hardware
You can find the hardware info in the file hw.info. Tests were conducted on a 6-core machine.




