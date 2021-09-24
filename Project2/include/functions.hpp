//declaring all the functions that will be used in project 2
//and including the headers needed



//include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>

//Creates a tridiagonal matrix
//n determines size of matrix
//l determines lower diagonal elements
//d determines diagonal elements
//u determines upper diagonal elements
arma::mat tridiag_matrix(int N, double l, double d, double u);

//calculates the analytical solution to eigenvalues and eigenvectors    
void analytical_solution(double a, double d, int N, arma::vec& eigval, arma::mat& eigvec);

//finds and saves the indicies of the largest off-diagonal element.
//Returns the value of that element
double max_offdiag_symmetric(const arma::mat& A, int &k, int &l);

// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);



#endif // end of include guard