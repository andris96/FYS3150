#include "functions.hpp"

int main()
{
    arma::mat A = {{1, 0, 0, 0.5}, {0, 1, -0.7, 0}, {0, -0.7, 1, 0}, {0.5, 0, 0, 1}};
    int k; 
    int l;

    std::cout << "Matrix A is given as:" << std::endl;
    A.print();
    std::cout << "The maximum off-diagonal value of A is: " << max_offdiag_symmetric(A,k,l) << std::endl;

    std::cout << "Which is the element with indicies k = " << k << ", l = " << l << std::endl; 



}