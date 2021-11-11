#ifndef __IsingModel_hpp__
#define __IsingModel_hpp__

#include <stdlib.h>
#include <string.h>
#include <armadillo>
#include <map>

class IsingModel {

    private:
    friend class TestIsingModel;
    int L;
    double T;
    arma::mat s;

    double kB;
    double beta;

    public:
    // Constructor
    IsingModel(int L_in, double T_in);

    // Generate a random spin configuration
    arma::mat generate_random_spin_config(int L);

    // Flip the spin at lattice site (i,j)
    void flip_spin(int i, int j);

    // Compute the energy difference due to a spin flip at site (i,j)
    int compute_energy_diff_due_to_flip(int i, int j);

    // Compute the map from deltaE to exp(-B*deltaE)
    std::map<int, double> make_energy_map();

    // Run a MC cycle using Markov Chains
    void run_MCMC(int N, std::string outputFileName);



};

#endif