#ifndef __TestIsingModel_hpp__
#define __TestIsingModel_hpp__

#include "IsingModel.hpp"
#include <armadillo>
#include <assert.h>

class TestIsingModel{
    private:
    //This class is a friend class of IsingModel
    //don't need anything that isn't already in IsingModel yet

    public:
    TestIsingModel(); // Empty constructor

    // Run all test
    void run_all_tests();

    // Calculating the different analytical solutions
    double Partition(double beta);
    double epsilon_squared_expectation(double beta, double N);
    double epsilon_expectation(double beta, double N);
    double m_squared_expectation(double beta, double N);
    double E_expectation(double beta);
    double E_squared_expectation(double beta);
    double M_squared_expectation(double beta);
    double M_abs_expectation(double beta);
    double Cv(double beta, double T, int N);
    double khi(double beta, int N);

    // Methods for testing member methods of class IsingModel
    void test_rand_uniform();
    void test_generate_random_spin_config();
    void test_initiate();
    void test_compute_energy_diff_due_to_flip();

};

#endif