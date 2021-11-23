#include <IsingModel.hpp>

// Constructor
IsingModel::IsingModel(int L_in, double T_in) {
    L = L_in;
    T = T_in;

    s = arma::Mat<int>(L, L);
    E = 0;
    M = 0;

    kB = 1.0; // Set proper value!!!
    beta = 1/(kB*T);

    energy_map = make_energy_map(); 
}

/**
 * Random number generator drawn from a uniform distrution using system clock as seed
 * 
 * Adapted from: https://github.com/anderkve/FYS3150/blob/master/code_examples/random_number_generation/main_minimal.cpp
 * (...)
 */
double IsingModel::rand_uniform() {
    unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator;
    generator.seed(seed);
    std::uniform_real_distribution<double> my_distribution(0.0, 1.0);
    return my_distribution(generator);
}

/**
 * Update the spin state s with a random configuration
 * 
 * The matrix s is initiated with all ones, and then each elem
 * if flipped with an equal probability. 
 * 
 * (...)
 * 
 */
void IsingModel::generate_random_spin_config() {
    s = arma::Mat<int>(L, L, arma::fill::ones);
    double r;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            r = rand_uniform();
            if (r > 0.50) {
                s(i, j) = -1;
            }
        }
    }
}

void IsingModel::generate_ordered_spin_config() {
    s = arma::Mat<int>(L, L, arma::fill::ones);
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
 * Note on the use of modulo operators
 * -----------------------------------
 * As we only are considering the neighbours below and to the right, we will not get any cases where
 * the modulo operator is used on negative numbers, hence no additional complications are needed as in 
 * compute_energy_diff_due_to_flip().
 * 
 */
void IsingModel::initiate(bool random) {
    if (random == true){
        generate_random_spin_config();
    }
    else{
        generate_ordered_spin_config();
    }
    M = arma::accu(s);
    E = 0;
    int current; int bottom; int right; // spin values
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            current = s(i, j);
            bottom = s((i + 1) % L, j);
            right = s(i, (j + 1) % L);
            E += (current*bottom + current*right);
        }
    }
    E *= -1;
}

// Set a spin config
void IsingModel::set_s(arma::Mat<int> s_in) {
    assert((s_in.n_rows == L) && (s_in.n_cols == L));
    s = s_in;
}


// Get the boltzmann constant
double IsingModel::get_Kb() {
    return kB;
}

// Calculate Cv
double IsingModel::Cv_calculate(double beta, int N, double T, double mean_E2, double mean_E){
    return beta/(N*T) * (mean_E2 - mean_E*mean_E);
}

// Calculate x
double IsingModel::X_calculate(double beta, int N, double T, double mean_M2, double mean_M_abs){
    return beta/N * (mean_M2 - mean_M_abs*mean_M_abs);
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
 * Note on the use of modulo operators
 * -----------------------------------
 * Due to how the % operator is implemented in C++, taking the modulo of a negative number
 * will yield a negative number, which is not what we want in order to include periodic 
 * boundary conditionds. So instead of taking a % b we need to take (b + (a%b)) % b. 
 * 
 * Source:
 * https://stackoverflow.com/questions/7594508/modulo-operator-with-negative-values
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
    int top = s((L + ((i - 1)%L)) % L, j); // In ordinary math notation: s(i-1 % L, j)
    int bottom = s((i + 1) % L, j);
    int left = s(i, (L + ((j - 1)%L)) % L); // s(i, j-1 % L)
    int right = s(i, (j + 1) % L);

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
 * max_trials (int) : The maximum number of trials to search for state with lower energy
 */
void IsingModel::metropolis(int max_trials) {
    int i; int j;
    int deltaE;
    for (int trial = 0; trial < max_trials; trial++) {

        // s.print();
        // std::cout << "\nE = " << E << "\n";

        // Pick random spin site from the lattice and flip it
        i = int(rand_uniform()*10*L) % L; //std::rand() % L;
        j = int(rand_uniform()*10*L) % L;//std::rand() % L;
        flip_spin(i, j);

        // Compute the energy difference due to the spin flip and ratio, but 
        // accept the new state immediately if deltaE <= 0
        deltaE = compute_energy_diff_due_to_flip(i, j);
        // if (deltaE <= 0) {
        //     accepted = true;
        //     break;
        // }
        double ratio = energy_map[deltaE];

        // Acceptance step: 
        // reject if r >= ratio, that is revert to previous state, ie. flip back the spin at (i,j)
        double r = rand_uniform();
        if (r >= ratio) {
            flip_spin(i, j);
        }
        else {
            // Update energy and magnetization'
            E += deltaE; //
            M += 2*s(i, j); // Magnetization changes by either +2 or -2
        }
    }
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
 * random (bool) : Initiates a random configuration if set to true, ordered config if false
 * samples (bool) : Creates a file called samples.txt to store each value of epsilon, if set to true
 */
void IsingModel::monte_carlo(int max_cycles, int max_trials, arma::vec &results, bool random, 
                             bool samples, const char* filename) {
    initiate(random);
    for (int cycle = 0; cycle < max_cycles; cycle++) {

        // Search for a lower energy/ higher probability state..
        metropolis(max_trials);

        // Store the relavant quantities
        results(0) += E;
        results(1) += E*E;
        results(2) += M;
        results(3) += M*M;
        results(4) += abs(M);

        // Saving epsilon in a file called samples.txt
        if (samples == true){
            int N = L*L;

            std::ofstream file;
            file.open(filename, std::ios_base::app);
            file << E/N << "\n";
            file.close();
        }
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
 * random (bool) : Initiate a random state if true, initiate an ordered state if false
 * print (bool) : if true, prints the various quantities to the terminal
 * expectation (bool) : If true, writes the expectation value of epsilon and m to text files
 * e_file (const char*) : Name of the file for <epsilon> values 
 * m_file (const char*) : Name of the file for <m> values
 */
void IsingModel::estimate_quantites_with_MCMC(int max_cycles, int max_trials,  arma::rowvec& evalues, 
                                              bool random, bool print, bool expectation,
                                              const char* e_file, const char* m_file) {
    arma::vec results = arma::vec(5, arma::fill::zeros);
    

    // Only interested in expectation values, therefore samples = false
    monte_carlo(max_cycles, max_trials, results, random, false);

    // Total number of spins
    int N = L*L;

    // Compute expectation values
    double mean_E = results(0)/max_cycles; 
    double mean_e = mean_E/N;
    double mean_E2 = results(1)/max_cycles;
    double mean_e2 = mean_E2/N;
    double mean_M = results(2)/max_cycles;
    double mean_m = mean_M/N;
    double mean_M2 = results(3)/max_cycles;
    double mean_m2 = mean_M2/N;
    double mean_M_abs = results(4)/max_cycles;
    double mean_m_abs = mean_M_abs/N;

    // Compute specific heat capacity Cv and magnetic susceptibility X per spin
    double Cv = Cv_calculate(beta, N, T, mean_E2, mean_E);
    double X = X_calculate(beta, N, T, mean_M2, mean_M_abs);

    evalues(0) = mean_e;
    evalues(1) = mean_m_abs;
    evalues(2) = Cv;
    evalues(3) = X;

    // Printing values to the terminal
    if (print == true){
        std::cout << "<e>: " << mean_e << std::endl;
        std::cout << "<e^2>: " << mean_e2 << std::endl;
        std::cout << "<|m|>: " << mean_m_abs << std::endl;
        std::cout << "<m^2>: " << mean_m2 << std::endl;
        std::cout << "Cv: " << Cv << std::endl;
        std::cout << "X: " << X << std::endl;
    }
    
    // Saving expectation values to files
    if (expectation == true){
        std::ofstream files;
    
        files.open(e_file, std::ios_base::app); // std::ios_base::app is to append the values
        files << mean_e << "\n";
        files.close();

        files.open(m_file, std::ios_base::app);
        files << mean_m_abs << "\n";
        files.close();
    }


}