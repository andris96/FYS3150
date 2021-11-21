#include <IsingModel.hpp>
#include <TestIsingModel.hpp>

#include "omp.h"
#include <iostream>
#include <fstream>

/**
 * To build:
 * g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo
 * 
 * To run:
 * ./main.exe
 */

/**
 * 
 * (detailed doc string..)
 * 
 * Program 0 (default) : Run tests for class IsingModel
 * Program 1 : Comparing numerical for L=2 against analytical (Problem 4)
 * ...
 * Program 4 : Estimate <e>, <|m|>  Cv and X for different L's and T's using OpenMP (Problem 8)
 * 
 */
int main() {

    std::cout << "Welcome to Project 4!\n\n"
    << "Choose one of the following programs to run:\n"
    << "0 : Run all test for the class IsingModel\n"
    << "1 : Compare the estimated quantities for a lattice of size L=2 against\n"
    << "the analytical solution\n"
    // ...
    << "4 : (Problem 8)\n\n"
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

        L2.estimate_quantites_with_MCMC(max_cycles, max_trials, random, true, false);

        // Analytical for 2x2
        std::cout << "\n<e> analytical: " << analytical.epsilon_expectation(beta,N) << std::endl;
        std::cout << "<e^2> analytical: " << analytical.epsilon_squared_expectation(beta,N) << std::endl;
        std::cout << "<m^2> analytical: " << analytical.m_squared_expectation(beta,N) << std::endl;
        std::cout << "Cv analytical: " << analytical.Cv(beta, T,  N) << std::endl;
        std::cout << "X analytical: " << analytical.khi(beta, N) << std::endl;
        break;
    }

    // Problem 8 : Estimatiation of quantities through parallellization with OpenMP 
    case 4: {
        for (int L = 40;  L <= 100; L += 20) {
            //...
            continue;
        }

    }

    default:
        std::cout << "No valid input found, terminating!" << std::endl;
        break;

    } // end switch

    return 0;
}

