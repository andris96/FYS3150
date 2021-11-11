#include <IsingModel.hpp>

// Constructor
IsingModel::IsingModel(int L_in, double T_in) {
    L = L_in;
    T = T_in;

    s = arma::mat(L, L);
    E = 0;
    M = 0;

    kB = 1.0; // Set proper value!!!
    beta = 1/(kB*T);

    energy_map = make_energy_map(); 
}

/**
 * Update the spin state s with a random configuration
 * 
 * (...)
 * 
 */
void IsingModel::generate_random_spin_config() {
    s = arma::mat(L, L, arma::fill::randu);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (s(i, j) > 0.5) {
                s(i, j) = 1;
            }
            else {
                s(i, j) = -1;
            }
        }
    }
}

/**
 * Initiate with a new spin configuration and compute the associated energy and magnetization
 * 
 * (...)
 * 
 * To avoid double counting of spin-pairs when computing the system energy, only the spin pairs
 * current<->right and current<->bottom are considered for site (i,j). Boundary conditions are taken
 * into consideration by the modulo operations.
 * 
 *  |c|  |r|   o    o
 *  
 *  |b|   o    o    o
 * 
 *   o    o    o    o
 * 
 *   o    o    o    o 
 * 
 */
void IsingModel::initiate() {
    generate_random_spin_config();
    M = arma::accu(s);
    E = 0;
    int current; int bottom; int right; // spin values
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            bottom = s(i + 1 % L, j);
            right = s(i, j + 1 % L);
            E += (current*bottom + current*right);
        }
    }
    E = -E;
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
 *     +  -  +  ->    +  +  +         Total system energy change = -4
 *        +              +
 * 
 * 3.     -        3'    -
 *     +  -  -  ->    +  +  -         Total system energy change = 0
 *        +              +
 * 
 * 4.     -        4'    -
 *     +  -  -  ->    +  +  -         Total system energy change = 4
 *        -              -
 * 
 * 5.     -        5'    -
 *     -  -  -  ->    -  +  -         Total system energy change = 8
 *        -              -
 * 
 * Hence, if after a single spin flip, the config. around (i,j) is eg. 5', then
 * we know that the Energy of all pairs around (i,j) has gone from E_before = 4x -1 = -4
 * to E_after = 4x +1 = 4, giving an energy difference of deltaE = +8. 
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
 */
std::map<int, double> IsingModel::make_energy_map() {
    std::map<int, double> energy_map = {
        {-8, exp(-beta*(-8))}, 
        {-4, exp(-beta*(-4))},
        {0, 1}, 
        {4, exp(-beta*(4))}, 
        {8, exp(-beta*(8))}
    };
    return energy_map;
}

/**
 * Metropolis acceptance steps: run a single Monte Carlo cycle
 * 
 * ... (Detailed docstring..)
 * 
 * Params
 * ------
 * max_trials (int) : The maxium number of trials to search for state with lower energy
 */
void IsingModel::metropolis(int max_trials, std::map<int, double> energy_map) {
    int i; int j;
    int deltaE;
    for (int trial = 0; trial < max_trials; trial++) {

        // Pick random spin site from the lattice and flip it
        i = std::rand() % L;
        j = std::rand() % L;
        flip_spin(i, j);

        // Compute the energy difference due to the spin flip and ratio, but 
        // accept the new state immediately if deltaE <= 0
        deltaE = compute_energy_diff_due_to_flip(i, j);
        if (deltaE <= 0) {
            break;
        }
        double ratio = energy_map[deltaE];

        // Acceptance step: 
        // reject if r >= ratio, that is revert to previous state, ie. flip back the spin at (i,j)
        double r = std::rand() % 1;
        if (r >= ratio) {
            flip_spin(i, j);
        }
    } // end for-loop

    // Update energy and magnetization
    E += deltaE;
    M += 2*s(i, j); // why x2? (copied the updating procedure in the compendium..)
}

/**
 * Monte Carlo computaion.. stores values for relavant quantities...
 * 
 * (Detailed docstring..)
 * 
 * Params
 * ------
 * max_cycle (int) : The maximum number of cycles to search for a stationary state..(?)
 * max_trials (int) : The maximum number of trials to search for a state with lower energy
 * results (arma::vec) : Vector for storing results, taken as a reference, shape (5, 1).
 *      The format is [E, E*E, M, M*M, |M|]. 
 */
void IsingModel::monte_carlo(int max_cycles, int max_trials, arma::vec &results) {
    initiate();
    for (int cycle = 0; cycle < max_cycles; cycle++) {

        // Search for a lower energy state..
        metropolis(max_trials, energy_map);

        // Store the relavant quantities
        results(0) += E;
        results(1) += E*E;
        results(2) += M;
        results(3) += M*M;
        results(4) += abs(M);
    }
}

/**
 * Estimate relavant quantities...
 * 
 * (...)
 * 
 * Params
 * ------
 * max_cycle (int) : The maximum number of cycles to search for a stationary state..(?)
 * max_trials (int) : The maximum number of trials to search for a state with lower energy..(?)
 */
void IsingModel::estimate_quantites_with_MCMC(int max_cycles, int max_trials) {
    arma::vec results = arma::vec(5, arma::fill::zeros);
    monte_carlo(max_cycles, max_trials, results);

    // Total number of spins
    int N = L*L;

    // Compute expectation values per spin
    double mean_e = results(0)/max_cycles/N;
    double mean_e2 = results(1)/max_cycles/N;
    double mean_m = results(2)/max_cycles/N;
    double mean_m2 = results(3)/max_cycles/N;
    double mean_m_abs = results(4)/max_cycles/N;

    // Compute specific heat capacity Cv and magnetic susceptibility X per spin
    double Cv = beta * (mean_e2 - mean_e*mean_e);
    double X = beta * (mean_m2 - mean_m*mean_m);

    // Do something more... print to terminal ... save to file for later plotting etc.. 
    // ...
    // ...
}