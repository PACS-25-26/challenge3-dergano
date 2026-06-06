#include <iostream>
#include <mpi.h>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>
#include "muParser.h"


double uexact(double x, double y);

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (argc != 5)
        {
        if (rank == 0)
            std::cerr << "Usage: " << argv[0] << " <matrix_size> <num_threads> <function> <boundary_function>" << std::endl;
        MPI_Finalize();
        return EXIT_FAILURE;
        }    
    
    int k;
    int num_threads;
    if (rank == 0){
        k = std::stoi(argv[1]);
        num_threads = std::stoi(argv[2]);
    }
    MPI_Bcast(&num_threads, 1, MPI_INT, 0, MPI_COMM_WORLD);
    omp_set_num_threads(num_threads);

    std::string f = argv[3];
    std::string boundary = argv[4];
    

    std::vector<std::size_t> num(5);
    for (int i = 0; i < 5; ++i){
        num[i] = std::pow(2, i + 4);
    }

    std::size_t max_iters = 100000;
    double tol = 1e-6;
    double error = 0.0;
    std::vector<double> errors;
    

    for (auto &n : num){

        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();    
        if (rank == 0){
            std::cout << "Running for n = " << n << std::endl;
        }
            
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double h = 1.0 / (n - 1);
        int rest = n % size;
        int local_rows = n / size + (rank < rest ? 1 : 0);

        std::size_t local_start = n * rank * (n / size) + n * std::min(rank, rest);

        std::vector<double> local_Uk(local_rows * n, 0.0);
        std::vector<double> local_Uk1(local_rows * n, 0.0);

        std::vector<double> prerow_Uk(n, 0.0);
        std::vector<double> postrow_Uk(n, 0.0);

        // Evaluation of forcing and boundary function, both given by the user (Parallel implementation)
        std::vector<double> local_f(local_rows * n, 0.0); 
        std::vector<double> local_boundary(local_rows * n, 0.0);
        std::vector<double> local_true_solution(local_rows * n, 0.0);
        #pragma omp parallel 
        {
            double thread_x = 0.0;
            double thread_y = 0.0;
            mu::Parser thread_parser_forcing;
            mu::Parser thread_parser_boundary;

            thread_parser_forcing.DefineVar("x", &thread_x);
            thread_parser_forcing.DefineVar("y", &thread_y);
            thread_parser_forcing.DefineConst("pi", M_PI);
            thread_parser_forcing.SetExpr(f);

            thread_parser_boundary.DefineVar("x", &thread_x);
            thread_parser_boundary.DefineVar("y", &thread_y);
            thread_parser_boundary.DefineConst("pi", M_PI);
            thread_parser_boundary.SetExpr(boundary);


            #pragma omp for
            for (std::size_t j = 0; j < local_rows * n; ++j) 
            {
                thread_x = (j % n) * h;
                thread_y = (local_start / n + j / n) * h;
                std::size_t col = j % n;
                std::size_t global_row = local_start / n + j / n;
                local_true_solution[j] = uexact(thread_x, thread_y);

                if (global_row == 0 || global_row == n - 1 || col == 0 || col == n - 1) {
                    local_boundary[j] = thread_parser_boundary.Eval();
                }
                else {
                    local_f[j] = thread_parser_forcing.Eval();
                }      
            }
        }

        double conv_crit = tol + 1.0;   
        double local_tol = 1e-4;
        int max_local_iters = 100;
        
        std::vector<double> prev_outer_Uk(local_rows * n, 0.0);

        for (std::size_t iter = 0; iter < max_iters; ++iter)
        {
            #pragma omp parallel for
            for (std::size_t j = 0; j < local_rows * n; ++j) {
                prev_outer_Uk[j] = local_Uk[j];
            }

            // 2. Scambio delle righe al bordo (Ghost Nodes fissi per il ciclo interno)
            if (rank > 0) {
                MPI_Sendrecv(local_Uk.data(), n, MPI_DOUBLE, rank - 1, 0,
                            prerow_Uk.data(), n, MPI_DOUBLE, rank - 1, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank < size - 1) {
                MPI_Sendrecv(local_Uk.data() + (local_rows - 1) * n, n, MPI_DOUBLE, rank + 1, 0,
                            postrow_Uk.data(), n, MPI_DOUBLE, rank + 1, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            //risolutore locale (Jacobi) con criteri di convergenza locale
            double local_inner_crit = local_tol + 1.0;
            int inner_iter = 0;

            while (local_inner_crit > local_tol && inner_iter < max_local_iters)
            {
                local_inner_crit = 0.0;

                #pragma omp parallel for 
                for (std::size_t j = 0; j < local_rows * n; ++j)
                {
                    std::size_t col = j % n;
                    std::size_t global_row = local_start / n + j / n;

                    if(global_row == 0 || global_row == n - 1 || col == 0 || col == n - 1) {
                        local_Uk1[j] = local_boundary[j]; 
                    }
                    else {
                        double down = (j < n) ? prerow_Uk[col] : local_Uk[j - n];
                        double up   = (j >= (local_rows - 1) * n) ? postrow_Uk[col] : local_Uk[j + n];
                        double left  = local_Uk[j - 1];
                        double right = local_Uk[j + 1];

                        local_Uk1[j] = 0.25 * (left + right + up + down + (h * h) * local_f[j]); 
                    }
                }
                
                // Errore locale interno (per capire la stabilità)
                #pragma omp parallel for reduction(+:local_inner_crit) 
                for (std::size_t j = 0; j < local_rows * n; ++j) {
                    local_inner_crit += std::pow(local_Uk1[j] - local_Uk[j], 2);
                } 
                local_inner_crit = std::sqrt(local_inner_crit);

                std::swap(local_Uk, local_Uk1);
                inner_iter++;
            }


            // controllo di convergenza
            double local_conv_crit = 0.0;
            double local_error = 0.0;

            // Confronto la soluzione attuale con prev_outer_Uk, non con la iterazione locale precedente!
            #pragma omp parallel for reduction(+:local_conv_crit, local_error) 
            for (std::size_t j = 0; j < local_rows * n; ++j)
            {
                local_conv_crit += std::pow(local_Uk[j] - prev_outer_Uk[j], 2);
                local_error += std::pow(local_Uk[j] - local_true_solution[j], 2);
            } 

            MPI_Allreduce(&local_conv_crit, &conv_crit, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&local_error, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            

            conv_crit = std::sqrt(h * conv_crit);
            error = std::sqrt(h * error);
            
            errors.push_back(error); 

            if (conv_crit < tol)
            {
                if (rank == 0) {
                    std::cout << "Global convergence reached in " << iter + 1 << " external iterations." << std::endl;
                }
                break; 
            }
        }

        if (rank == 0 && conv_crit >= tol)
        {
            std::cout << "Did not converge in " << max_iters << " iterations. Final convergence criterion: " << conv_crit << std::endl;
        }

        std::vector<double> global_Uk;
        std::vector<int> recvcounts(size);
        std::vector<int> displs(size);

        if (rank == 0)
        {
            global_Uk.resize(n * n); 
            int current_j = 0;
            
            for (int i = 0; i < size; ++i)
            {
                int rows_for_rank_i = (n / size) + (i < rest ? 1 : 0);
                recvcounts[i] = rows_for_rank_i * n;
                displs[i] = current_j;
                current_j += recvcounts[i];
            }
        }

        MPI_Gatherv(local_Uk.data(), local_rows * n, MPI_DOUBLE,
                    global_Uk.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        if (rank == 0)
        {
            std::ofstream out("laplace_solution.vtk");
            
            //Intestazione obbligatoria
            out << "# vtk DataFile Version 3.0\n";
            out << "Laplace equation solution\n";
            out << "ASCII\n";
            
            //Definizione della griglia spaziale
            out << "DATASET STRUCTURED_POINTS\n";
            out << "DIMENSIONS " << n << " " << n << " 1\n"; 
            out << "ORIGIN 0 0 0\n";
            out << "SPACING " << h << " " << h << " 1\n";
            
            // I dati della soluzione
            out << "POINT_DATA " << n * n << "\n";
            out << "SCALARS u double 1\n";
            out << "LOOKUP_TABLE default\n";
            
            // Scrittura del vettore globale
            for (std::size_t i = 0; i < n * n; ++i)
            {
                out << global_Uk[i] << "\n";
            }
            
            out.close();
            std::cout << "Risultato esportato in laplace_solution.vtk" << std::endl;
        }  
        if (rank == 0) {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            std::cout << "Parallel execution completed for n = " << n << std::endl;
            std::cout << "Elapsed time: " << elapsed_seconds.count() << "s\n";
        }          
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}


double uexact(double x, double y)
{
    return std::sin(2 * M_PI * x) * std::sin(2 * M_PI * y);
}