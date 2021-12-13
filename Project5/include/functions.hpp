#include <armadillo>
#include <iostream>
#include <iomanip> 
#include <complex>
#include "math.h"

using namespace std::complex_literals;  

int convertk(int i, int j, int M);

void AB(int M, double h, double dt, arma::mat V, arma::sp_cx_mat &A, arma::sp_cx_mat &B);

arma::cx_vec u_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py);

arma::cx_vec u_inner_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py);

arma::cx_mat vec_to_mat(arma::cx_vec u);

arma::mat V_config(int slits, double v0, arma::vec x, arma::vec y);

