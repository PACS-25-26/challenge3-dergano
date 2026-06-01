#include<iostream>
#include "matrix.hpp"
#include <mpi.h>
#include <chrono>
#include <string>
#include <vector>
#include <cmath>


double f(double x, double y);

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
    
        int n;
    if (rank == 0){
        std::cout << "Matrix size: " << argv[1] << std::endl;
        n = std::stoi(argv[1]);
        std::vector<double> Uk(n * n, 0.0);
        std::vector<double> Uk1(n * n, 0.0);

        // Initialize Uk with initial conditions
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                Uk[i * n + j] = 0.0; // Initial condition, can be modified in the extra
            }
        }
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double h = 1.0 / (n - 1);
    int rest = n % size;
    int local_rows = n / size + (rank < rest ? 1 : 0);
    std::size_t local_start = rank * (n / size) + n * std::min(rank, rest);
    std::size_t local_end = local_start + local_rows * n;
    std::vector<double> local_Uk(local_rows * n, 0.0);
    std::vector<double> local_Uk1(local_rows * n, 0.0);
    std::vector<double> prerow_Uk(n, 0.0);
    std::vector<double> postrow_Uk(n, 0.0);

    // Scatter the initial Uk to all processes
    MPI_Scatterv(rank == 0 ? Uk.data() : nullptr, 
                 std::vector<int>(size, local_rows * n).data(), 
                 std::vector<int>(size, 0).data(), 
                 MPI_DOUBLE, 
                 local_Uk.data(), 
                 local_rows * n, 
                 MPI_DOUBLE, 
                 0, 
                 MPI_COMM_WORLD);

    int max_iters = 10000;
    double tol = 1e-6;
    double err = tol+1;
    double sum = 0.0;

    for (std::size_t iter = 0; iter < max_iters; ++iter)
    {
        // Scambio delle righe al bordo
        if (rank > 0)
        {
            MPI_Send(local_Uk(local_Uk.size() - n).data(), 
                    n, 
                    MPI_DOUBLE, 
                    rank + 1, 
                    0, 
                    MPI_COMM_WORLD);
            MPI_Recv(prerow_Uk.data(),
                     n, 
                     MPI_DOUBLE, 
                     rank - 1, 
                     0, 
                     MPI_COMM_WORLD, 
                     MPI_STATUS_IGNORE);
        }
        if (rank < size - 1)
        {
            MPI_Send(local_Uk.data(), 
                    n, 
                    MPI_DOUBLE, 
                    rank - 1, 
                    0, 
                    MPI_COMM_WORLD);
            MPI_Recv(postrow_Uk.data(),
                     n, 
                     MPI_DOUBLE, 
                     rank + 1, 
                     0, 
                     MPI_COMM_WORLD, 
                     MPI_STATUS_IGNORE);
        }

        // Calcolo di Uk1

        for (std::size_t j = 0; j < local_row * n; ++j)
        {
            double x = (j % n) * h;
            double y = (local_start / n + j / n) * h;
            if (local_start < n || j % n == 0 || j % n == n - 1 || local_end > (n - 1) * n)
            {
                for (auto it = Uk1.begin(); it != Uk1.end(); ++it)
                    Uk1(it)= Uk(it); // Boundary condition, can be modified in the extra
            }
            else 
            {
                if (j < n)
                {
                Uk1(j) = (Uk(j-1) + Uk(j+1) + prerow_Uk(j) + Uk(j+n) + (h * h)*f(x, y));
                }
                else if (j >= (local_rows - 1) * n)
                {
                    Uk1(j) = (Uk(j-1) + Uk(j+1) + Uk(j-n) + postrow_Uk(j) + (h * h)*f(x, y));
                }
                else 
                {
                    Uk1(j) = (Uk(j-1) + Uk(j+1) + Uk(j-n) + Uk(j+n) + (h * h)*f(x, y));
                }   
            }
        }
        // Calcolo dell'errore
        double local_sum = 0.0;
        for (std::size_t j = 0; j < local_rows * n; ++j)
        {
            local_sum += std::pow(Uk1(j) - Uk(j), 2);
        }

        if (rank == 0)
        {
            MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            err = std::sqrt(h * sum);
            if (err < tol)
            {
                std::cout << "Converged in " << iter + 1 << " iterations with error: " << err << std::endl;
                break;
            }
        }
        

    }


        




    

    
    MPI_Finalize();
    return EXIT_SUCCESS;
}


double f(double x, double y)
{
    return 8 * (M_PI * M_PI) * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
}
