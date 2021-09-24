#include "functions.hpp"
#include <cmath>

//creates a tridiagonal matrix
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


//calculating the largest off-diagonal element in the symmetrical matrix A
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


// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

    double tau = (A(l,l) - A(k,k))/(2*A(k,l));
    //tan cos and sin theta
    double t;
    double c;
    double s;

    if(tau > 0){
        t = 1/(tau + sqrt(1 + std::pow(tau,2)));
    }
    else{
        t = 1/(tau - sqrt(1 + std::pow(tau,2)));
    }

    c = 1/(sqrt(1 + std::pow(tau,2)));
    s = c*t;

    //transforming matrix elements of A
    double old_akk = A(k,k); //need to save A^{m}_{kk} to not calculate A_{ll} with A^{m+1}_{kk}
    //updating these 4 elements explicitly as they differ from the general rule
    A(k,k) = A(k,k)*std::pow(c,2) - 2*A(k,l) + A(l,l)*std::pow(s,2);
    A(l,l) = A(l,l)*std::pow(c,2) + 2*A(k,l) + old_akk*std::pow(s,2);
    A(k,l) = 0;
    A(l,k) = 0;

    //declaring variables for the upcoming for-loop
    int N = A.n_rows;
    double old_aik; //This is needed for the same reason as old_aik above

    //Updating the rest of the elements in A
    for(int i = 0; i<N; i++){
        //iterating through every element where i != k,l as they have already been set above
        if((i != k) && (i != l)){
            old_aik = A(i,k);
            A(i,k) = A(i,k)*c - A(i,l)*s;
            A(k,i) = A(i,k);
            A(i,l) = A(i,l)*c + old_aik;
            A(l,i) = A(i,l);
        }
    }

    //updating elements in R
    for(int i = 0; i<N; i++){
        double old_rik = R(i,k);
        R(i,k) = R(i,k) - R(i,l);
        R(i,l) = R(i,l) + old_rik;
    }
}

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    int k;
    int l;
    int N = A.n_cols;
    arma::mat A_rotated = A;
    eigenvectors.set_size(N,N);
    eigenvectors = arma::mat(N,N,arma::fill::eye); //R^{1} = I
    eigenvalues.set_size(N);
    
    iterations = 0;


    while(max_offdiag_symmetric(A_rotated,k,l) < eps){
        max_offdiag_symmetric(A_rotated,k,l);
        //std::cout << k << std::endl;
        jacobi_rotate(A_rotated,eigenvectors,k,l);
        iterations++;
        //std::cout << iterations << std::endl;
        if(iterations == maxiter){
            std::cout << "maximum number of iterations was reached" << std::endl;
            break;
        converged = true;
        }
    }

    A_rotated.print();
    /*
    for(int i; i<N; i++){
        eigenvalues(i) = A_rotated(i,i);
    }
    */
    //arma::sort(eigenvectors);
    //arma::sort(eigenvalues);

}