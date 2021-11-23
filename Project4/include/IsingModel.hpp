#ifndef __IsingModel_hpp__
#define __IsingModel_hpp__

#include <stdlib.h>
#include <string.h>
#include <random>
#include <chrono>
#include <assert.h>
#include <armadillo>
#include <map>

/**
 * Class for ...
 * 
 * (Detailed docstring..)
 * 
 * (The structure is partly based on the program/script at p. 436-437 in the compendium/ document p. 448-449)
 * 
 * 
 */
class IsingModel {

    private:
    friend class TestIsingModel;
    int L;
    double T;

    arma::Mat<int> s;
    double E;
    double M;

    double kB;
    double beta;

    std::map<int, double> energy_map;

    public:
    // Constructor
    IsingModel(int L_in, double T_in);

    // Random number generators with system clock as seed
    double rand_uniform();

    // Update the spin state with a random configuration
    void generate_random_spin_config();

    // Update the spin state with an ordered configuration
    void generate_ordered_spin_config();

    // Setters (only for testing purposes)
    void set_s(arma::Mat<int> s_in); 


    // Getters
    double get_Kb();

    // Calculate Cv and x
    double Cv_calculate(double beta, int N, double T, double mean_E2, double mean_E);
    double X_calculate(double beta, int N, double T, double mean_M2, double mean_M_abs);

    // Initiate with a new spin configuration and compute the associated energy and magnetization
    // If random == true, the spin configuration is randomized
    // If random == false, the spin configuration is ordered
    void initiate(bool random);

    // Flip the spin at lattice site (i,j)
    void flip_spin(int i, int j);

    // Compute the energy difference due to a spin flip at site (i,j)
    int compute_energy_diff_due_to_flip(int i, int j);

    // Compute the map from deltaE to exp(-B*deltaE)
    std::map<int, double> make_energy_map();

    // Metropolis/ acceptance step: find a new state with lower E
    void metropolis(int max_trials);

    // Monte Carlo computations...
    // If random == true, it generates a random spin configuration
    // If random == false, it generates an ordered configuration
    void monte_carlo(int max_cycles, int max_trials, arma::vec &results, bool random, bool samples,
                     const char* filename = "samples.txt");

    // Estimate relavant quantities...
    void estimate_quantites_with_MCMC(int max_cycles, int max_trials, arma::rowvec evalues, 
                                      bool random = false, bool print = false, bool expectation = false, 
                                      const char* e_file = "e_file.txt", 
                                      const char* m_file = "m_file.txt");

};

#endif