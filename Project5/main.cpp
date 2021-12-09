#include <iostream>
#include <iomanip>      
#include <armadillo>
#include <complex>
#include <cmath>

using namespace std::complex_literals;


/**
 * 
 * To build:
 * g++ main.cpp -o main.exe -larmadillo
 * 
 * To run:
 * ./main.exe
*/


// Temporary setup, gonna move functions over to a src file at some point

/**
 * Takes indices i and j of a (M-2)x(M-2) matrix and convert them to a corresponding vector index k
*/
int convertk(int i, int j, int M){
  //We let K start from 1 and not 0, since j and i start at 1 aswell (We are excluding boundary conditions)
  return ((M-2)*(j-1) + i);
}



/**
 * Constructs matrices A, and assign the values of vector a as diagonal in A
 * Also takes value r and fills the matrices according to the Crank-Nicolson scheme in two dimensions
 * To construct matrix B, all you have to do is input another vector b and -r instead of a and r
*/
arma::cx_mat A(std::complex<double> r, arma::cx_vec a){

  int L = a.size();

  arma::cx_mat A, B;
  A.set_size(L,L);

  A(0,0) = a(0);
  A(0,1) = -r;
  A(L-1,L-2) = -r;
  A(L-1,L-1) = a(L-1); 

  // Setting up a tridiagonal matrix
  for(int i=1; i<L-1; i++){
        A(i, i) = a(i);
        A(i, i-1) = -r; 
        A(i, i+1) = -r;
  }

  // Setting up the rest of the diagonal elements
  int s = sqrt(L);
  for(int i = 0; i < (L-s); i++){
    A(i,s+i) = -r;
    A(s+i,i) = -r;

  }

  arma::cx_mat matrix;
  matrix = A;
  return matrix;
  
}

//Not tested yet
arma::cx_mat construct_matrix(int M, double h, double dt, arma::mat V){

  arma::cx_vec a(M-2);
  arma::cx_vec b(M-2);

  int length = (M-2)*(M-2);
  int k;
  std::complex<double> r= 1i*dt/(2*h*h);
  
  for(int i = 0; i < length; i++){
    for(int j = 0; j < length; j++){
      k = convertk(i,j,M);
      a(k) = 1. + 4.*r + 1i * dt/2. * V(i,j);
      b(k) = 1. - 4.*r - 1i * dt/2. * V(i,j);
    }
  }

  return A(r,a);
}



int main()
{

  int i = 4;
  int j = 3;
  int M = 7;

  int K = convertk(i,j,M);
  std::cout << K << std::endl;

  arma::cx_vec a(9, arma::fill::ones);
  arma::cx_vec b = a + a;
  double r = 4;

  arma::cx_mat mat = A(r,a);
  arma::mat A = arma::real(mat);

  std::cout << std::setw(3) << A << std::endl;
 


  return 0;
}