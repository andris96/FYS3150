#include "functions.hpp"
#include <cmath>

//to run: g++ problem3_main.cpp src/functions.cpp 
//-I include -o problem3_main.exe -larmadillo -llapack -lblas


int main()
{
    int n = 7;
    int N = n-1;
    double stepsize = 1/double(n);



    double a = -1/std::pow(stepsize,2.0);
    double b = 2/std::pow(stepsize,2.0);

    arma::mat A = tridiag_matrix(N,a,b,a);
    
    arma::vec eigval;
    arma::mat eigvec;
    arma::vec analytical_eigval;
    arma::mat analytical_eigvec;
    


    arma::eig_sym(eigval,eigvec,A);

    analytical_solution(a,b,N,analytical_eigval,analytical_eigvec);

    arma::mat norm_eigvec = arma::normalise(eigvec);
    arma::mat norm_analytical_eigvec = arma::normalise(analytical_eigvec);

    std::cout << "eigvalues:" << std::endl;
    eigval.print();
    //analytical_eigval.print();
    //norm_eigvec.print();
    //norm_analytical_eigvec.print();


    
}