#include <IsingModel.hpp>


/**
 * To build:
 * g++ burn_in.cpp src/IsingModel.cpp -I include -o burn_in.exe -larmadillo
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
    files.open("T1R_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T2O_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T2R_e_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T1O_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T1R_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T2O_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("T2R_m_values.txt", std::ofstream::out | std::ofstream::trunc);
    files.open("cycles.txt", std::ofstream::out | std::ofstream::trunc);
    files.close();

    arma::ivec max_cycles = arma::regspace<arma::ivec>(100, 100, 1000);
    int max_trials = 1000;

    files.open("cycles.txt");
    files << max_cycles << std::endl;
    files.close();

    bool random = true;
    bool ordered = false;

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