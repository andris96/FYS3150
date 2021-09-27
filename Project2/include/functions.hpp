//declaring all the functions that will be used in project 2
//and including the headers needed



//include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <armadillo>

//creates a tridiagonal matrix with l as lower diagonal
//d as diagonal, and u as upperdiagonal elements.
arma::mat tridiag_matrix(int N, double l, double d, double u);

//Function calculating the analytical eigenvalues and vectors of a symmetric tridiagonal matrix
//with d as diagonal and a as lower and upper-diagonal  
void analytical_solution(double a, double d, int N, arma::vec& eigval, arma::mat& eigvec);

//finds and saves the indicies of the largest off-diagonal element.
//Returns the value of that element, and stores the indecies in k and l
double max_offdiag_symmetric(const arma::mat& A, int &k, int &l);

//Function for doing a single jacobi rotation
//assumes that the input matrix is symmetric
//Takes the matrix A, which is to be solved
//Also takes a rotation matrix R as input
//k and l are the indicies for the largest off-diagonal element
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

//This function uses the function jacobi_rotate to rotate matrix A multiple times.
//The largest off-diagonal element is found by utilizing max_offdiag_symmetric.
//It takes in matrix A, epsilon which is the limit for how big off-diagonal elements are allowed to be
//Also takes a vector to place the eigenvalues, and a matrix to place eigenvectors
//It takes in iterations, which should start at 0, and will change the value for each iteration
//maxiter is the maximum number of iterations we allow 
//converged is a variable that will tell us wether or not the maximum number of iterations was reached 
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);



#endif // end of include guard