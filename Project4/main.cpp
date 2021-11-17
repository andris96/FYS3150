#include <IsingModel.hpp>
#include <TestIsingModel.hpp>

// To build: 
// g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo


int main() {

    // tmp: run all tests
    TestIsingModel tests = TestIsingModel();
    tests.run_all_tests();

    // Set params
    double T = 1.0;
    int L = 2;
    double beta = 1./T;

    int max_trials = 10000;
    int max_cycles = 500;

    // Instantiate
    IsingModel L2(L, T);
    TestIsingModel analytical;

    L2.estimate_quantites_with_MCMC(max_cycles, max_trials);

    // Analytical for 2x2
    std::cout << "Cv analytical: " << analytical.Cv(beta, T,  max_trials) << std::endl;
    std::cout << "X analytical: " << analytical.khi(beta, max_trials) << std::endl;

    return 0;
}

