#include "functions.hpp"
#include <iostream>

int main(){
    int n = 7;
    int N = n-1;
    double stepsize = 1/double(n);
    double a = -1/std::pow(stepsize,2.0);
    double b = 2/std::pow(stepsize,2.0);
    arma::mat A = tridiag_matrix(N,a,b,a);
    arma::mat R = arma::mat(N,N,arma::fill::eye);

    /*
    A.print();

    jacobi_rotate(A,R,0,1);

    A.print();
    */
    
    std::cout << "potato" << std::endl;
    double eps = 0.0001;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1000;
    int iterations;
    bool converged = false;

    jacobi_eigensolver(A,eps, eigenvalues, eigenvectors, maxiter,iterations, converged);

    eigenvectors.print();
    eigenvalues.print();
    
}
