#include <IsingModel.hpp>
#include <TestIsingModel.hpp>

/**
 * To build:
 * g++ main.cpp src/IsingModel.cpp src/TestIsingModel.cpp -I include -o main.exe -larmadillo
 * 
 * To run:
 * ./main.exe
 */


int main() {

    // tmp: run all tests
    TestIsingModel tests = TestIsingModel();
    tests.run_all_tests();

    // Set params
    double T = 1.6;
    int L = 2;
    int N = L*L;
    double beta = 1./T;
    bool random = true;

    // User input
    int max_trials, max_cycles;
    std::cout << "Max trials: ";
    std::cin >> max_trials;
    std::cout << "Max cycles: ";
    std::cin >> max_cycles;
    std::cout << "\n";

    // Instantiate
    IsingModel L2(L, T);
    TestIsingModel analytical;

    L2.estimate_quantites_with_MCMC(max_cycles, max_trials, random);

    // Analytical for 2x2
    std::cout << "\n<e> analytical: " << analytical.epsilon_expectation(beta,N) << std::endl;
    std::cout << "<e^2> analytical: " << analytical.epsilon_squared_expectation(beta,N) << std::endl;
    std::cout << "<m^2> analytical: " << analytical.m_squared_expectation(beta,N) << std::endl;
    std::cout << "Cv analytical: " << analytical.Cv(beta, T,  N) << std::endl;
    std::cout << "X analytical: " << analytical.khi(beta, N) << std::endl;


    return 0;
}

