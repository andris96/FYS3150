#include "functions.hpp"
#include <cmath>

//to compile: g++ problem3_main.cpp src/functions.cpp -I include -o problem3_main.exe -larmadillo
//to run: ./problem3_main.exe

int main()
{
    int n = 7; //amount of steps
    int N = n-1; //Size of matrix
    double stepsize = 1/double(n);

    //Elements of the tridiagonal matrix are given by a and d
    //a is the upper- and lower-diagonal elements, d are the diagonal elements
    double a = -1/std::pow(stepsize,2.0);
    double d = 2/std::pow(stepsize,2.0);
    arma::mat A = tridiag_matrix(N,a,d,a);
    
    //Declaring matrices and vectors for eigenvalues and eigenvectors
    arma::vec eigval;
    arma::mat eigvec;
    arma::vec analytical_eigval;
    arma::mat analytical_eigvec;
    //using eig_sym
    arma::eig_sym(eigval,eigvec,A);

    //calculating the analytical solution
    analytical_solution(a,d,N,analytical_eigval,analytical_eigvec);

    //normalizing the eigenvectors to compare them
    arma::mat norm_eigvec = arma::normalise(eigvec);
    arma::mat norm_analytical_eigvec = arma::normalise(analytical_eigvec);

    std::cout << "Matrix A: \n" << A << std::endl;

    std::cout << "Eigenvalues from eig_sym: \n";
    eigval.print();
    std::cout << "Analytical eigenvalues: \n";
    analytical_eigval.print();

    std::cout << "Normalised eigenvectors from eig_sym: \n";
    norm_eigvec.print();
    std::cout << "Normalised analytical eigenvectors: \n";
    norm_analytical_eigvec.print();


    
}