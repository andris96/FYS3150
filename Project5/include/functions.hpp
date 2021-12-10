#include <armadillo>
#include <iostream>
#include <iomanip> 
#include <complex>

using namespace std::complex_literals;  

int convertk(int i, int j, int M);

void AB(int M, std::complex<double> h, std::complex<double> dt, arma::cx_mat V, arma::cx_mat &A, arma::cx_mat &B);

