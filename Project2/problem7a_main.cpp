#include "functions.hpp"
#include <iostream>
#include <iomanip>

//to compile: g++ problem7a_main.cpp src/functions.cpp -I include -o problem7a_main.exe -larmadillo
//to run: ./problem7a_main.exe
//Alot of repeating code in the main programs, but it at least it works :) 
int main(){
    int n = 10; 
    int N = n-1;
    double stepsize = 1/double(n);
    double a = -1/std::pow(stepsize,2.0);
    double b = 2/std::pow(stepsize,2.0);
    arma::mat A = tridiag_matrix(N,a,b,a);
    
    double eps = 1e-7;
    arma::vec eigenvalues(N, arma::fill::zeros);
    arma::mat eigenvectors = arma::mat(N,N,arma::fill::eye);
    int maxiter = 1e6;
    int iterations = 0;
    bool converged = true; 

    jacobi_eigensolver(A,eps, eigenvalues, eigenvectors, maxiter,iterations, converged);

    if(converged == true){
        std::cout << "convergence was reached before hitting maximum number of iterations." << std::endl;
    }
    else{
        std::cout << "Maximum number of iterations was reached." << std::endl;
    }

    std::ofstream myfile; 
    myfile.open ("eigenvec10.txt"); 
    myfile << "eigenvectors\n";

    for (int i = 0; i<N; i++ ){
    myfile << std::scientific << std::setprecision(5) << eigenvectors(i,0) << " " << eigenvectors(i,1)
    << " " << eigenvectors(i,2) << "\n";
    }
    

    
}