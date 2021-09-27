#include "functions.hpp"
#include <iostream>

//to compile: g++ problem6_main.cpp src/functions.cpp -I include -o problem6_main.exe -larmadillo
//to run: ./problem6_main.exe

//this program could be improved and more dynamic, but it does what I need it to.
//I could make the output be a file, which then can be read by a python program to plot.
//I could also make a better function for making the tridiagonal matrices 
//such that I wouldn't have to use a for-loop with so many steps.
 
int main(){
    //defining some values which are independent of N
    double eps = 1e-7;
    int maxiter = 1e+4;

    //Using a for-loop to do the same as problem5_main.cpp with varying sizes of N
    for (int i = 10; i < 71; i += 5){
        //Defining N, a and b
        int N = i;
        double stepsize = 1/double(N+1);
        double a = -1/std::pow(stepsize,2.0);
        double b = 2/std::pow(stepsize,2.0);

        //making the matrices needed
        arma::mat A = tridiag_matrix(N,a,b,a);
        arma::vec eigenvalues(N, arma::fill::zeros);
        arma::mat eigenvectors = arma::mat(N,N,arma::fill::eye);

        //start on 0 iterations, and converged is set to true
        int iterations = 0;
        bool converged = true; 

        jacobi_eigensolver(A,eps, eigenvalues, eigenvectors, maxiter,iterations, converged);

        //Printing the results
        std::cout << "Number of iterations: " << iterations << std::endl;

    }

    
    
}