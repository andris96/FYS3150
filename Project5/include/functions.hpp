#include <armadillo>
#include <iostream>
#include <iomanip> 
#include <complex>

using namespace std::complex_literals;  

int convertk(int i, int j, int M);

void AB(int M, std::complex<double> h, std::complex<double> dt, arma::cx_mat V, arma::cx_mat &A, arma::cx_mat &B);

std::complex<double> u_init(std::complex<double> x, std::complex<double> y, std::complex<double> xc, std::complex<double> yc, 
                            std::complex<double> sx, std::complex<double> sy, std::complex<double> px, std::complex<double> py);

arma::cx_mat V_config(int slits, std::complex<double> v0, std::complex<double> slitsize);
