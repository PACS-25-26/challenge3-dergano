### RESULT DISCUSSION
We are discussing the result reported in the scalability_report.txt file. (To generate it you may rely on the instruction in the README file)

## TABLE OF RESULTS
| N | Serial (1 P, 1 T) | MPI (2 P, 1 T) | MPI (4 P, 1 T) | OpenMP (1 P, 2 T) | OpenMP (1 P, 4 T) | Hybrid (2 P, 2 T) |  Hybrid (3 P, 2 T) 

### Time
|---------|--------|--------|--------|---------|---------|--------|---------|
| **16**  | 0.0009 | 0.0018 | 0.0016 | 0.0061  | 0.0166  | 0.0256 | 0.3120  |
| **32**  | 0.0050 | 0.0592 | 0.0067 | 0.0237  | 0.0619  | 0.0615 | 0.2583  |
| **64**  | 0.0489 | 0.0772 | 0.0461 | 0.1199  | 0.2640  | 0.1029 | 0.8239  |
| **128** | 0.6847 | 0.3850 | 0.4386 | 0.9668  | 1.5173  | 0.6793 | 2.5512  |
| **256** | 9.8055 | 5.3011 | 4.1569 | 11.0254 | 13.1570 | 6.3954 | 16.7058 |
*(P = MPI Processes / Ranks, T = OpenMP Threads)*

### Iterations
|---------|------|------|------|------|------|------|------|
| **16**  | 37   | 21   |   64 | 37   | 37   | 21   | 69   |
| **32**  | 144  | 26   |  214 | 144  | 144  | 26   | 207  |
| **64**  | 524  | 376  | 731  | 524  | 524  | 376  | 692  |
| **128** | 1838 | 2118 | 2520 | 1838 | 1838 | 2118 | 2367 |
| **256** | 6240 | 7369 | 8698 | 6240 | 6240 | 7369 | 8080 |
*(P = MPI Processes / Ranks, T = OpenMP Threads)*

### Speed-up (S)
| N | MPI (2 P, 1 T) | MPI (4 P, 1 T) | OpenMP (1 P, 2 T) | OpenMP (1 P, 4 T) | Hybrid (2 P, 2 T) | Hybrid (3 P, 2 T) |
|---------|------|------|------|------|------|------|
| **16**  | 0.50 | 0.56 | 0.15 | 0.05 | 0.04 |~0.00 |
| **32**  | 0.08 | 0.75 | 0.21 | 0.08 | 0.08 | 0.02 |
| **64**  | 0.63 | 1.06 | 0.41 | 0.19 | 0.48 | 0.06 |
| **128** | 1.78 | 1.56 | 0.71 | 0.45 | 1.01 | 0.27 |
| **256** | 1.85 | 2.36 | 0.89 | 0.75 | 1.53 | 0.59 |
*(P = MPI Processes / Ranks, T = OpenMP Threads)*

### Efficiency (E)
| N | MPI (2 P, 1 T) | MPI (4 P, 1 T) | OpenMP (1 P, 2 T) | OpenMP (1 P, 4 T) | Hybrid (2 P, 2 T) | Hybrid (3 P, 2 T) |
|---------|------|------|------|------|------|------|
| **16**  | 0.25 | 0.14 | 0.08 | 0.01 | 0.01 |~0.00 |
| **32**  | 0.04 | 0.19 | 0.11 | 0.02 | 0.02 |~0.00 |
| **64**  | 0.32 | 0.27 | 0.21 | 0.05 | 0.12 | 0.01 |
| **128** | 0.89 | 0.39 | 0.36 | 0.11 | 0.25 | 0.05 |
| **256** | 0.93 | 0.59 | 0.45 | 0.19 | 0.38 | 0.10 |
*(P = MPI Processes / Ranks, T = OpenMP Threads)*


## COMMENTS
We observed that the fastest approach to solving the problem is a pure MPI parallelization. This result was expected; since the algorithm involves relatively simple mathematical operations, it is more efficient to execute them serially at the local level rather than splitting the workload across multiple OpenMP threads. In fact, the overhead introduced by managing and synchronizing the threads outweighs the computational benefits.

The number of external iterations also aligns with theoretical expectations. Using OpenMP parallelization does not alter the total number of iterations, as it merely parallelizes local computations without changing the underlying numerical method. Conversely, we observe a significant variation in iterations when using MPI, which is a direct consequence of the domain decomposition approach (Block Jacobi).

Finally, the speed-up and efficiency metrics confirm our observations, with one notable detail: solving the problem with 2 MPI processes yields a higher parallel efficiency compared to using 4 MPI processes.