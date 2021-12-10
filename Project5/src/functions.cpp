#include "functions.hpp"

/**
 * Takes indices i and j of a (M-2)x(M-2) matrix and convert them to a corresponding vector index k
*/
int convertk(int i, int j, int M){
    //We let K start from 1 and not 0, since j and i start at 1 aswell (We are excluding boundary conditions)
    return ((M-2)*(j-1) + i);
}





/**
 * Constructs matrices A and B according to the Crank-Nicolson scheme in two dimensions
 * 
*/
void AB(int M, double h, double dt, arma::mat V, arma::cx_mat &A, arma::cx_mat &B){

    int L = (M-2)*(M-2);
    std::complex<double> r = 1.i*dt/(2*h*h);

    arma::cx_vec a(L, arma::fill::zeros);
    arma::cx_vec b(L, arma::fill::zeros);

    int k = 0;

    //Calculating the elements of a and b
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            k = convertk(i,j,M);
            a(k) = 1. + 4.*r + 1.i * dt/2. * V(i,j);
            b(k) = 1. - 4.*r - 1.i * dt/2. * V(i,j);
            }
    }

    A(0,0) = a(0);
    A(0,1) = -r;
    A(L-1,L-2) = -r;
    A(L-1,L-1) = a(L-1); 

    B(0,0) = b(0);
    B(0,1) = r;
    B(L-1,L-2) = r;
    B(L-1,L-1) = b(L-1); 


    // Setting up a tridiagonal matrix
    for(int i=1; i<L-1; i++){
        A(i, i) = a(i);
        A(i, i-1) = -r; 
        A(i, i+1) = -r;

        B(i, i) = b(i);
        B(i, i-1) = r; 
        B(i, i+1) = r;
    }


        
    // Setting up the rest of the off-diagonal elements
    int s = sqrt(L);
    for(int i = 0; i < (L-s); i++){
        A(i,s+i) = -r;
        A(s+i,i) = -r;

        B(i,s+i) = -r;
        B(s+i,i) = -r;

    }
  
}