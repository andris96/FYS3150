#include<armadillo>

// Solve the problem for one particle..
arma::mat solve_analytical_1p(double v0, double x0, double z0, double tmax, int steps, int q, 
                              double B0, double V0, double m, double d);