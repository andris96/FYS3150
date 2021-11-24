#include<burn_in.hpp>

/**
 * (Doc string..)
 */
void estimate_burn_in() {
    
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

    // // For storing expectation values of each max_cycles
    // // row format: { max_cycles T1Oe T1Om T1Re T1Rm T2Oe T2Om T2Re T2Rm}
    // arma::mat expectation_values(max_cycles_step_size, 9, arma::fill::zeros);

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
}