#include <IsingModel.hpp>
#include <TestIsingModel.hpp>

#include "omp.h"
#include <iostream>
#include <fstream>

/**
 * In order to run in parallel, add -fopenmp flag when building the program
 * 
 * To build:
 * g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo
 * 
 * To run:
 * ./main.exe
 */

/**
 * 
 * 
 * Program 0 (default) : Run tests for class IsingModel
 * Program 1 : Comparing numerical for L=2 against analytical (Problem 4)
 * Program 2 : Make estimates of <epsilon> and <|m|> under different circumstances
 * Program 3 : Estimate probability distributions
 * Program 4 : Estimate <e>, <|m|>  Cv and X for different L's and T's using OpenMP (Problem 8)
 * 
 */
int main() {

    std::cout << "Welcome to Project 4!\n\n"
    << "Choose one of the following programs to run:\n"
    << "0 : Run all tests for the class IsingModel\n"
    << "1 : Compare the estimated quantities for a lattice of size L=2 against\n"
    << "    the analytical solution\n"
    << "2 : Make estimates of <epsilon> and <|m|> under different circumstances\n" 
    << "    in order to estimate burn-in time\n"
    << "3 : Generate samples of epsilon to estimate the probability distributions\n"
    << "4 : Estimate <e>, <|m|>  Cv and X for different L's and T's\n\n"
    << "These programs generate text files containing values. Plotting of these \n"
    << "values are done in a separate python program\n\n"
    << "Enter an integer between 0-4: ";

    int program;
    std::cin >> program;
    switch (program) {

    // Run all tests
    case 0: {
        TestIsingModel tests = TestIsingModel();
        tests.run_all_tests();
        break;
    }

    // Compare numerical against anlytical for L=2 (Problem 4)
    // (COULD INCLUDE SOME PLOTTING IN ORDER TO ANSWER PROBLEM 4B, IF TIME)
    case 1: {
        // Set params
        double T = 1.6;
        int L = 2;
        int N = L*L;
        double beta = 1./T;
        bool random = true;

        // User input
        int max_trials, max_cycles;
        std::cout << "Max trials: ";
        std::cin >> max_trials;
        std::cout << "Max cycles: ";
        std::cin >> max_cycles;
        std::cout << "\n";

        // Instantiate
        IsingModel L2(L, T);
        TestIsingModel analytical;
        arma::rowvec expectation_values = arma::rowvec(4, arma::fill::zeros);
        
        L2.estimate_quantites_with_MCMC(max_cycles, max_trials, expectation_values, random, true, false);

        // Analytical for 2x2
        std::cout << "\n<e> analytical: " << analytical.epsilon_expectation(beta,N) << std::endl;
        std::cout << "<e^2> analytical: " << analytical.epsilon_squared_expectation(beta,N) << std::endl;
        std::cout << "<m^2> analytical: " << analytical.m_squared_expectation(beta,N) << std::endl;
        std::cout << "Cv analytical: " << analytical.Cv(beta, T,  N) << std::endl;
        std::cout << "X analytical: " << analytical.khi(beta, N) << std::endl;
        break;
    }

    // Finding <epsilon> and <|m|> with different number of MC cycles,
    // In order to estimate burn-in time. (Problem 5)
    // Plotting is done in a separate python program
    case 2: {

        // Set parameters
        double T1 = 1.0;
        double T2 = 2.4; 
        int L = 20;
        int max_trials = 1000;

        int max_cycles_start = 3000;
        int max_cycles_end = 5000;
        int max_cycles_step_size = 100;

        // Parmeters for estimate_quantites_with_MCMC()
        bool random = true;
        bool ordered = false;
        arma::rowvec temp(4,arma::fill::zeros);

        // Start timer
        auto start = std::chrono::steady_clock::now();
        
        // Parallelization (build with -fopenmp)
        #ifdef _OPENMP
        {
        #pragma omp parallel
        {
        // Each thread gets unique output file
        int my_thread = omp_get_thread_num();
        std::ofstream ofile;
        std::string filename = "expectation_values_thread_" + std::to_string(my_thread) + ".txt";
        ofile.open(filename.c_str(), std::ofstream::trunc);

        arma::rowvec temp(4,arma::fill::zeros);

        #pragma omp for
        for (int max_cycles = max_cycles_start; max_cycles <= max_cycles_end; max_cycles += max_cycles_step_size){
            IsingModel T1_ordered(L,T1); 
            IsingModel T1_random(L,T1); 
            IsingModel T2_ordered(L,T2);
            IsingModel T2_random(L,T2);

            ofile << max_cycles << " ";

            T1_ordered.estimate_quantites_with_MCMC(max_cycles, max_trials, temp, ordered, false, false);
            ofile << temp(0) << " " << temp(1) << " "; //T1Oe T1Om

            T1_random.estimate_quantites_with_MCMC(max_cycles, max_trials, temp, random, false);
            ofile << temp(0) << " " << temp(1) << " "; //T1Re T1Rm

            T2_ordered.estimate_quantites_with_MCMC(max_cycles, max_trials, temp, ordered, false);
            ofile << temp(0) << " " << temp(1) << " "; //T2Oe T1Om

            T2_random.estimate_quantites_with_MCMC(max_cycles, max_trials, temp, random, false);
            ofile << temp(0) << " " << temp(1) << std::endl; //T2Re T1Rm
        }
        ofile.close();
        } // end omp parallel
        } // end ifdef

        #else
        {
        IsingModel T1_ordered(L,T1); 
        IsingModel T1_random(L,T1); 
        IsingModel T2_ordered(L,T2);
        IsingModel T2_random(L,T2);

        arma::rowvec e = arma::rowvec(5, arma::fill::zeros);
        for (int max_cycles = max_cycles_start; max_cycles <= max_cycles_end; max_cycles += max_cycles_step_size){
            T1_ordered.estimate_quantites_with_MCMC(max_cycles, max_trials, e, ordered, false,
                                                    true,"T1O_e_values.txt","T1O_m_values.txt" );
            T1_random.estimate_quantites_with_MCMC(max_cycles, max_trials, e, random, false,
                                                true, "T1R_e_values.txt","T1R_m_values.txt");
            T2_ordered.estimate_quantites_with_MCMC(max_cycles, max_trials, e, ordered, false,
                                                    true, "T2O_e_values.txt","T2O_m_values.txt");
            T2_random.estimate_quantites_with_MCMC(max_cycles, max_trials, e, random, false,
                                                true, "T2R_e_values.txt","T2R_m_values.txt");
            }
        }
        #endif

        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Time spent on estimating: " << elapsed_seconds.count() << "s\n";
        break;
    }

    // Generating samples of epsilon at two different temperatures with ordered
    // and random states (Problem 6). 
    // Plotting is done in a separate python program. 
    case 3: {

        auto start = std::chrono::steady_clock::now();

        // Set parameters (same as case 2)
        double T1 = 1.0;
        double T2 = 2.4; 
        int L = 20;
        int max_trials = 1000;

        IsingModel T1samples(L,T1);
        IsingModel T2samples(L,T2);
        int max_cycles = 5000;
        arma::vec resultsT1 = arma::vec(5, arma::fill::zeros);
        arma::vec resultsT2 = arma::vec(5, arma::fill::zeros);

        std::ofstream files;
        files.open("samplesT1.txt", std::ofstream::trunc);
        files.close();
        files.open("samplesT2.txt", std::ofstream::trunc);
        files.close();

        // We could parallelize this code to do one thread for each temperature, however this seemed
        // to not save much time, and was troublesome to implement as we don't always get 
        // the amount of threads we request for.
        T1samples.monte_carlo(max_cycles, max_trials, resultsT1, true, true, "samplesT1.txt");
        T2samples.monte_carlo(max_cycles, max_trials, resultsT2, true, true, "samplesT2.txt");
        
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Program run-time: " << elapsed_seconds.count() << "s\n";

        break;
    }

    // Problem 8 : Estimatiation of quantities through parallellization with OpenMP 
    case 4: {
        int max_cycles = 1000;
        int max_trials = 1000;
        arma::vec T = arma::linspace(2.3, 2.4, 10);
        arma::mat expectation_values(T.size(), 4, arma::fill::zeros);
        arma::rowvec temp(4,arma::fill::zeros);
        
        auto start = std::chrono::steady_clock::now();
        for (int L = 40;  L <= 100; L += 20) { 
            #ifdef _OPENMP
            {
            #pragma omp parallel for
            for (int i = 0; i < T.size(); i++){
                IsingModel LT(L, T(i)); 
                LT.estimate_quantites_with_MCMC(max_cycles, max_trials, temp);
                expectation_values.row(i) = temp;
                        
            }
            expectation_values.save("expectation_valuesL" + std::to_string(L) + ".txt", arma::raw_ascii);
            }
            #else
            {
            for (int i = 0; i < T.size(); i++){
                IsingModel LT(L, T(i)); 
                LT.estimate_quantites_with_MCMC(max_cycles, max_trials, temp);
                expectation_values.row(i) = temp;
            }

            expectation_values.save("expectation_valuesL" + std::to_string(L) + ".txt", arma::raw_ascii);
            }
            #endif
            
        }
        
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Time spent on estimating: " << elapsed_seconds.count() << "s\n";

        break;
    }

    default:
        std::cout << "No valid input found, terminating!" << std::endl;
        break;

    } // end switch


    return 0;
}

