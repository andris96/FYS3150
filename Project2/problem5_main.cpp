#include "functions.hpp"
#include <iostream>

//to compile: g++ problem5_main.cpp src/functions.cpp -I include -o problem5_main.exe -larmadillo
//to run: ./problem5_main.exe

int main(){
    //making matrix A once again
    int n = 7;
    int N = n-1;
    double stepsize = 1/double(n);
    double a = -1/std::pow(stepsize,2.0);
    double b = 2/std::pow(stepsize,2.0);
    arma::mat A = tridiag_matrix(N,a,b,a);
    
    //defining the parameters needed for jacobi_eigensolver
    double eps = 1e-7;
    arma::vec eigenvalues(N, arma::fill::zeros);
    arma::mat eigenvectors = arma::mat(N,N,arma::fill::eye);
    int maxiter = 1000;
    int iterations = 0;
    bool converged = true; 

    jacobi_eigensolver(A,eps, eigenvalues, eigenvectors, maxiter,iterations, converged);

    //making the program state wether or not the maximum number of iterations was reached
    if(converged == true){
        std::cout << "convergence was reached before hitting maximum number of iterations." << std::endl;
    }
    else{
        std::cout << "Maximum number of iterations was reached." << std::endl;
    }

    //printing iterations, eigenvectors and eigenvalues
    std::cout << "Number of iterations: " << iterations << std::endl;

    std::cout << "eigenvectors: " << std::endl;
    eigenvectors.print();
    std::cout << "eigenvalues: " << std::endl;
    eigenvalues.print();

    
}
