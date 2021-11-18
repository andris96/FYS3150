#include <IsingModel.hpp>

int main() {
    
    double T1 = 1.0;
    double T2 = 2.4; 
    int L = 20;
    int N = L*L;
    double beta1 = 1./T1;
    double beta2 = 1./T2;

    arma::ivec max_cycles = arma::regspace<arma::ivec>(0, 100, 2000);
    int max_trials = 1000;



    // Instantiate
    IsingModel T1_ordered(L,T1); 
    IsingModel T1_random(L,T1); 
    IsingModel T2_ordered(L,T2);
    IsingModel T2_random(L,T2);

    for (int i = 0; i < max_cycles.size(); i++){
        T1_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, false);
        T1_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, true);
        T2_ordered.estimate_quantites_with_MCMC(max_cycles(i), max_trials, false);
        T2_random.estimate_quantites_with_MCMC(max_cycles(i), max_trials, true);
    }
}