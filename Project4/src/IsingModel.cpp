#include <IsingModel.hpp>

// Constructor
IsingModel::IsingModel(int L_in, double T_in) {
    L = L_in;
    T = T_in;
    s = generate_random_spin_config(L);

    kB = 1.0; // Set proper value!!!
    beta = 1/(kB*T);
}

/**
 * Generate a random (lattice) spin configuration
 * 
 * Returns
 * -------
 * s_new (arma::mat) : A randomly generated spin configuration, shape (L, L)
 */
arma::mat IsingModel::generate_random_spin_config(int L) {
    arma::mat s_new = arma::mat(L, L, arma::fill::randu);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (s_new(i, j) > 0.5) {
                s_new(i, j) = 1;
            }
            else {
                s_new(i, j) = -1;
            }
        }
    }

}

// Flip the spin at lattice site (i,j)
void IsingModel::flip_spin(int i, int j) {
    s(i, j) *= -1;
}

/** 
 * Compute the energy difference due to a single spin flip at lattice site (i,j)
 * 
 * After a single spin flip, the configuration in the neighbourhood of (i,j)
 * will be in one of five (unique) spin configurations.
 * 
 * 1.     +        1'    +
 *     +  -  +  ->    +  +  +         Total system energy change = -8
 *        +              +
 * 
 * 2.     -        2'    -
 *     +  -  +  ->    +  +  +         Total system energy change = 
 *        +              +
 * 
 * 3.     -        3'    -
 *     +  -  -  ->    +  +  -         Total system energy change = 
 *        +              +
 * 
 * 4.     -        4'    -
 *     +  -  -  ->    +  +  -         Total system energy change = 
 *        -              -
 * 
 * 5.     -        5'    -
 *     -  -  -  ->    -  +  -         Total system energy change = 
 *        -              -
 * 
 * Hence, if after a single spin flip, the config. around (i,j) is eg. 5', then
 * we know that the Energy of all pairs around (i,j) has increased by 1, ie. that 
 * the Energy of the whole lattice has increased by 4x1 = 4. 
 * 
 * So knowing the local spin configuration around site (i,j) after a single spin flip,
 * tells us how much the system energy has changed from the previous state/ lattice spin
 * config.
 * 
 * ...
 * 
 * params
 * ------
 * i, j (int) : Lattice site (i,j)
 * 
 * returns
 * -------
 * (int) : The energy difference due to a single spin flip at site (i,j)
 * 
 */
int IsingModel::compute_energy_diff_due_to_flip(int i, int j) {
    // The spin values at site (i,j) before and after flip
    int site_before = -s(i, j);
    int site_after = s(i, j);

    // Extract spin values for the neighbourhood, with periodic BC taken into consideration
    int top = s(i - 1 % L, j);
    int bottom = s(i + 1 % L, j);
    int left = s(i, j - 1 % L);
    int right = s(i, j + 1 % L);

    // Compute energy contributions from the local neighbourhooud and return the difference
    int energy_before = -(site_before*top + site_before*bottom + site_before*left + site_before*right);
    int energy_after = -(site_after*top + site_after*bottom + site_after*left + site_after*right);
    return energy_after - energy_before;
}

/**
 * Make a mapping from deltaE -> exp(-B*deltaB)
 * 
 * Make a map for storing the exp. values associated with each possible deltaE..
 * 
 * NOT COMPLETE!!! FILL IN REST OF THE MAPPINGS..
 * 
 */
std::map<int, double> IsingModel::make_energy_map() {
    std::map<int, double> energy_map = {
        {-8, exp(-beta*(-8))}
    };
}

/**
 * Run a single Monte Carlo cycle
 * 
 * ... (Detailed docstring..)
 * 
 * Params
 * ------
 * N (int) : The number of repitations
 * outputFileName (std::string) : Name of the .txt-file in which the results are stored
 * 
 */
void IsingModel::run_MCMC(int N, std::string outputFileName) {

    // Create the energy map in order to avoid computing exp(.) in the loops
    std::map<int, double> energy_map = make_energy_map();

    int i;
    int j;
    arma::mat results = arma::mat(N, 4);
    for (int n = 0; n < N; n++) {

        // Pick random spin from lattice and flip it
        i = std::rand() % L;
        j = std::rand() % L;
        flip_spin(i, j);

        // Compute the energy difference due to the spin flip and 
        // the ratio p(s_after)/(p_s_before)
        int deltaE = compute_energy_diff_due_to_flip(i, j);
        double ratio = energy_map[deltaE];

        // Acceptance step: 
        // reject if r >= ratio, that is revert to previous state, ie. flip back the spin at (i,j)
        double r = std::rand() % 1;
        if (r >= ratio) {
            flip_spin(i, j);
        }

        // Compute relavant quantities...
        // ...
        // ...

        // Store the results in the matrix $results
        // ... 
        // ...
        // ...
    }
}