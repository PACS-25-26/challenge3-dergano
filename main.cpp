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
    }
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double h = 1.0 / (n - 1);
    int rest = n % size;
    int local_rows = n / size + (rank < rest ? 1 : 0);

    std::size_t local_start = n * rank * (n / size) + n * std::min(rank, rest);
    std::size_t local_end = local_start + local_rows * n;

    std::vector<double> local_Uk(local_rows * n, 0.0);
    std::vector<double> local_Uk1(local_rows * n, 0.0);

    std::vector<double> prerow_Uk(n, 0.0);
    std::vector<double> postrow_Uk(n, 0.0);

    std::size_t max_iters = 10000;
    double tol = 1e-6;
    double err = tol + 1.0;
    double sum = 0.0;



    for (std::size_t iter = 0; iter < max_iters; ++iter)
    {
        // Scambio delle righe al bordo (Ghost Nodes)
        if (rank > 0)
        {
            // Manda la prima riga al vicino di sopra, e ricevi la sua ultima riga dentro prerow
            MPI_Sendrecv(local_Uk.data(), n, MPI_DOUBLE, rank - 1, 0,
                         prerow_Uk.data(), n, MPI_DOUBLE, rank - 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1)
        {
            // Manda l'ultima riga al vicino di sotto, e ricevi la sua prima riga dentro postrow
            MPI_Sendrecv(local_Uk.data() + (local_rows - 1) * n, n, MPI_DOUBLE, rank + 1, 0,
                         postrow_Uk.data(), n, MPI_DOUBLE, rank + 1, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Calcolo di Uk1

        for (std::size_t j = 0; j < local_rows * n; ++j)
        {
            double x = (j % n) * h;
            double y = (local_start / n + j / n) * h;
            std::size_t col = j % n;

            if(local_start < n || col == 0 || col == n - 1 || local_end > (n - 1) * n)
            {
                local_Uk1[j]= 0.0; // Boundary condition, can be modified in the extra
            }
            else 
            {
                double down = (j < n) ? prerow_Uk[col] : local_Uk[j - n];

                // Il vicino "SOPRA" (y maggiore) si trova dopo in memoria (j + n) o nella postrow
                double up   = (j >= (local_rows - 1) * n) ? postrow_Uk[col] : local_Uk[j + n];

                double left  = local_Uk[j - 1];
                double right = local_Uk[j + 1];

                local_Uk1[j] = 0.25 * (left + right + up + down + (h * h) * f(x, y)); 
            }
        }
        // Calcolo dell'errore
        double local_sum = 0.0;
        for (std::size_t j = 0; j < local_rows * n; ++j)
        {
            local_sum += std::pow(local_Uk1[j] - local_Uk[j], 2);
            local_Uk[j] = local_Uk1[j]; // Piccolo trucco: aggiorna Uk già qui per risparmiare un intero ciclo for!
        }

        MPI_Allreduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        err = std::sqrt(h * sum);

        if (err < tol)
        {
            if (rank == 0) 
            {
                std::cout << "Converged in " << iter + 1 << " iterations with error: " << err << std::endl;
            }
            break; 
        }

        for (std::size_t j = 0; j < local_rows * n; ++j)
        {
            local_Uk[j]= local_Uk1[j];
        }

    }
    if (rank == 0 && err >= tol)
    {
        std::cout << "Did not converge in " << max_iters << " iterations. Final error: " << err << std::endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}


double f(double x, double y)
{
    return 8 * (M_PI * M_PI) * std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
}
