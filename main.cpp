#include<iostream>
#include "matrix.hpp"
#include <mpi.h>
#include <chrono>
#include <string>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 2)
        {
        if (rank == 0)
            std::cerr << "Usage: " << argv[0] << " <matrix_size>" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
        }
    
    unsigned int matrix_size = std::stoul(argv[1]);
    matrix       A(matrix_size, matrix_size);
    matrix       B(matrix_size, matrix_size);
    matrix       C(matrix_size, matrix_size);
    
    A.randomize();
    B.randomize();
    
    auto start = std::chrono::high_resolution_clock::now();
    C          = A * B;
    auto end   = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> elapsed = end - start;
    
    if (rank == 0)
        std::cout << "Time taken for matrix multiplication: " << elapsed.count() << " seconds" << std::endl;
    
    MPI_Finalize();
    return EXIT_SUCCESS;
}