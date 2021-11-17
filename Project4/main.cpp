#include <IsingModel.hpp>
#include <TestIsingModel.hpp>

// To build: 
// g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo


int main() {

    // tmp: run all tests
    TestIsingModel tests = TestIsingModel();
    tests.run_all_tests();

    // // Set params
    // double T = 1.0;
    // int L = 2;

    // int max_trials = 100;
    // int max_cycles = 20;

    // // Instantiate
    // IsingModel L2 = IsingModel(L, T);
    // L2.estimate_quantites_with_MCMC(max_cycles, max_trials);

    return 0;
}

