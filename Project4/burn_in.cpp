#include <IsingModel.hpp>
#include "omp.h"


/**
 * To build:
 * g++ burn_in.cpp src/IsingModel.cpp -I include -o burn_in.exe -larmadillo 
 * (add -fopenmp to run in parallel)
 * 
 * To run:
 * ./burn_in.exe
 */

int main() {
    
    double T1 = 1.0;
    double T2 = 2.4; 
    int L = 20;
    int N = L*L;
    double beta1 = 1./T1;
    double beta2 = 1./T2;


    // Instantiate
    IsingModel T1_ordered(L,T1); 
    IsingModel T1_random(L,T1); 
    IsingModel T2_ordered(L,T2);
    IsingModel T2_random(L,T2);
    
    // Creating empty files to save data
    std::ofstream files;
    files.open("T1O_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T1R_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T2O_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T2R_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T1O_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T1R_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T2O_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("T2R_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("cycles.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();

    arma::ivec max_cycles = arma::regspace<arma::ivec>(50, 50, 100);
    int max_trials = 1000;

    files.open("cycles.txt");
    files << max_cycles << std::endl;
    files.close();

    bool random = true;
    bool ordered = false;
    // In order to run in parallel, add -fopenmp flag when running the program
    #ifdef _OPENMP
    {
    #pragma omp parallel for
    for (int i = 0; i < max_cycles.size(); i++){
        T1_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, ordered, false,
                                                true,"T1O_e_values.txt","T1O_m_values.txt" );
        T1_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, random, false,
                                               true, "T1R_e_values.txt","T1R_m_values.txt");
        T2_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, ordered, false,
                                                true, "T2O_e_values.txt","T2O_m_values.txt");
        T2_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, random, false,
                                               true, "T2R_e_values.txt","T2R_m_values.txt");
        }
    }
    #else
    {
    for (int i = 0; i < max_cycles.size(); i++){
        T1_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, ordered, false,
                                                true,"T1O_e_values.txt","T1O_m_values.txt" );
        T1_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, random, false,
                                               true, "T1R_e_values.txt","T1R_m_values.txt");
        T2_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, ordered, false,
                                                true, "T2O_e_values.txt","T2O_m_values.txt");
        T2_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, random, false,
                                               true, "T2R_e_values.txt","T2R_m_values.txt");
        }
    }
    #endif
    // Doing problem 6 here, can be moved to main or somewhere else, or we can rename this file
    /*
    IsingModel T1samples(L,T1);
    IsingModel T2samples(L,T2);
    int cycles = 500;
    arma::vec expectation_values = arma::vec(5, arma::fill::zeros);

    files.open("samplesT1.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();
    files.open("samplesT2.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();

    T1samples.monte_carlo(cycles, max_trials, expectation_values, random, true, "samplesT1.txt");
    T2samples.monte_carlo(cycles, max_trials, expectation_values, random, true, "samplesT2.txt");
    */

    return 0;
}