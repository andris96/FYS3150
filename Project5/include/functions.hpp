#include <armadillo>
#include <iostream>
#include <iomanip> 
#include <complex>
#include "math.h"

using namespace std::complex_literals;  

/**
 * Takes indices i and j of a M x M matrix and convert them to a corresponding vector index k
*/
int convertk(int i, int j, int M);


/**
 * Constructs matrices A and B according to the Crank-Nicolson scheme in two dimensions
 * 
 * M is the number of points on x and y axis
 * s is the number of points on x and y axis excluding boundaries, so it is M-2
 * L is the number of rows and columns in A and B, which is the same as (M-2)*(M-2)
 * dt is the size of the timesteps
 * V is the potential matrix
 * A and B are the matrices we store the results in
 * a and b are the vectors containing the diagonal elements of A and B (not to be confused with the b = B*u in problem 3 and in main.cpp)
*/
void AB(int M, double h, double dt, arma::mat V, arma::sp_cx_mat &A, arma::sp_cx_mat &B);

/**
 * Constructs the initial wave function as a gaussian wave packet
 * Takes x and y coordinates and calculates the wave function in that point with given parameters
 * xc and yc are the center coordinates of the packet 
 * sx and sy are the widths of the packets in x and y direction
 * px and py are the momentums in x and y direction
*/
arma::cx_vec u_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py);

/** 
 * In order to solve the matrix equations with A and B, we make an almost identical function where the boundaries are excluded 
 * Takes x and y coordinates and calculates the wave function in that point with given parameters
 * xc and yc are the center coordinates of the packet 
 * sx and sy are the widths of the packets in x and y direction
 * px and py are the momentums in x and y direction
 *
 * This could be done in a different way by just altering the u_init to do both cases
*/
arma::cx_vec u_inner_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py);

/**
 * Takes a M^2-dimentional vector and converts it to a M x M matrix by utilizing the function convertk
 */
arma::cx_mat vec_to_mat(arma::cx_vec u);

/**
 * Constructs the V matrix at x and y coordinates with either 1 2 or 3 slits
 * v0 is the potential of the barrier
 * The size of each slit is set to 0.05
 * And the thickness is set to 0.02
 */
arma::mat V_config(int slits, double v0, arma::vec x, arma::vec y);

