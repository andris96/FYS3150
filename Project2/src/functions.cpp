#include "functions.hpp"
#include <cmath>

//creates a tridiagonal matrix with l as lower diagonal
//d as diagonal, and u as upperdiagonal elements.
arma::mat tridiag_matrix(int N, double l, double d, double u){

    //creating an identity matrix A of size NxN,
    arma::mat A = arma::mat(N,N,arma::fill::eye);

    //the first and last rows are exceptions to the general case
    //so they are defined before the for-loop
    A(0,0) = d;
    A(0,1) = u;
    A(N-1,N-2) = l;
    A(N-1,N-1) = d; 

    //start at the second row, end at the second to last row
    for(int i=1; i<N-1; i++){
        A(i,i) = d;
        A(i,i-1) = l;
        A(i,i+1) = u;
    }
    
    return A;
}

//Function calculating the analytical eigenvalues and vectors of a symmetric tridiagonal matrix
//with d as diagonal and a as lower and upper-diagonal
void analytical_solution(double a, double d, int N, arma::vec& eigval, arma::mat& eigvec){
    double pi = 2*acos(0.0); //calculating pi
    eigval.set_size(N); //setting the size of the vector and matrix
    eigvec.set_size(N,N);
    //I let i go from 0 to N-1 to get to all elements.
    //The calculation is based on the given formula, but the indicies are different
    //since the first element in the vector is i=0 and not i=1.
    for(int i = 0; i<N; i++){
        eigval(i) = d + 2*a*cos((i+1)*pi/(N+1)); 
        
        for(int j = 0; j<N; j++){
            eigvec(i,j) = sin((j+1)*(i+1)*pi/(N+1));
    
        }
        
    }

}

//finds and saves the indicies of the largest off-diagonal element.
//Returns the value of that element, and stores the indecies in k and l
double max_offdiag_symmetric(const arma::mat& A, int &k, int &l){
    int N = A.n_rows;
    k = 0;
    l = 1;

    double max_val = A(k,l);
    //Going through all columns
    for (int i = 0; i<N; i++){
        int j = i+1; //only checking the elements after A(i,i)
        while(j<N){
            if (std::abs(A(i,j)) > std::abs(max_val)){
                max_val = A(i,j);
                k = i;
                l = j;
            }
            j++;
        }
    }

    return max_val;

}


//Function for doing a single jacobi rotation
//assumes that the input matrix is symmetric
//Takes the matrix A, which is to be solved
//Also takes a rotation matrix R as input
//k and l are the indicies for the largest off-diagonal element
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

    double tau = (A(l,l) - A(k,k))/(2*A(k,l));
    //tan cos and sin theta
    double t;
    double c;
    double s;

    //defining tau
    if(tau > 0){
        t = 1/(tau + sqrt(1 + std::pow(tau,2)));
    }
    else{
        t = -1/(-tau + sqrt(1 + std::pow(tau,2)));
    }

    c = 1/(sqrt(1 + std::pow(t,2)));
    s = c*t;

    //transforming matrix elements of A
    double old_akk = A(k,k); //need to save A^{m}_{kk} to not calculate A_{ll} with A^{m+1}_{kk}
    //updating these 4 elements explicitly as they differ from the general rule
    A(k,k) = A(k,k)*std::pow(c,2) - 2*A(k,l)*c*s + A(l,l)*std::pow(s,2);
    A(l,l) = A(l,l)*std::pow(c,2) + 2*A(k,l)*c*s + old_akk*std::pow(s,2);
    A(k,l) = 0;
    A(l,k) = 0;

    //declaring variables for the upcoming for-loop
    int N = A.n_rows;
    double old_aik; //This is needed for the same reason as old_akk above

    //Updating the rest of the elements in A
    for(int i = 0; i<N; i++){
        //iterating through every element where i != k,l as they have already been set above
        if((i != k) && (i != l)){
            old_aik = A(i,k);
            A(i,k) = A(i,k)*c - A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + old_aik*s;
            A(l,i) = A(i,l);
        }
    }

    //updating elements in R
    for(int i = 0; i<N; i++){
        double old_rik = R(i,k);
        R(i,k) = R(i,k)*c - R(i,l)*s;
        R(i,l) = R(i,l)*c + old_rik*s;
    }
}

//This function uses the function jacobi_rotate to rotate matrix A multiple times.
//The largest off-diagonal element is found by utilizing max_offdiag_symmetric.
//It takes in matrix A, epsilon which is the limit for how big off-diagonal elements are allowed to be
//Also takes a vector to place the eigenvalues, and a matrix to place eigenvectors
//It takes in iterations, which should start at 0, and will change the value for each iteration
//maxiter is the maximum number of iterations we allow 
//converged is a variable that will tell us wether or not the maximum number of iterations was reached 
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    //k and l are defined in the function max_offdiag_symmetric
    int k; 
    int l;
    //Making the matrices/vectors which will contain eigenvectors and eigenvalues
    int N = A.n_cols;
    arma::mat A_rotated = A;
    eigenvectors.set_size(N,N);
    eigenvectors = arma::mat(N,N,arma::fill::eye); //R^{1} = I
    eigenvalues.set_size(N);
    
    //start at 0 iterations
    iterations = 0;
    while(std::abs(max_offdiag_symmetric(A_rotated,k,l)) > eps){
        //rotating once
        jacobi_rotate(A_rotated,eigenvectors,k,l);
        //increasing the rotation count by one
        iterations++;

        //Breaking the while loop after a certain number of iterations
        if(iterations == maxiter){
            converged = false;
            break;
        }
    }
     
    //The eigenvalues are set to the diagonal elements of A after the rotations
    for(int i = 0; i<N; i++){
        eigenvalues(i) = A_rotated(i,i);
    }

    //Ordering the eigenvectors and eigenvalues
    //only need to look at eigenvalues and put them in ascending order
    //then the corresponding column vectors can follow the same ordering
    arma::uvec indicies = arma::sort_index(eigenvalues);
    //need to store the old matrix to have a referance
    arma::vec eigenvalues_wrong_order = eigenvalues;
    arma::mat eigenvectors_wrong_order = eigenvectors;

    //changing the order of eigenvalues and eigenvectors which are column vectors
    for(int i = 0; i < N; i++){
        eigenvalues(i) = eigenvalues_wrong_order(indicies(i));
        eigenvectors.col(i) = eigenvectors_wrong_order.col(indicies(i));
    }
    

}