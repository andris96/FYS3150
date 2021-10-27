#include<armadillo>
#include "utils.hpp"

// Solve the problem for one particle..
arma::mat solve_analytical_1p(double v0, double x0, double z0, double tmax, int steps, int q, 
                              double B0, double V0, double m, double d) {
    arma::mat motion_r = arma::mat(steps, 3, arma::fill::zeros);

    double x = x0;
    double y = 0.0;
    double z = z0;
    double t = 0.0;
   
    double dt = tmax/steps;

    double omega_0 = q*B0/m;
    double omega_z_squared = 2*q*V0/(m*pow(d,2));
    double omega_z = sqrt(omega_z_squared);

    double omega_pluss = 1/2*(omega_0 + sqrt(pow(omega_0,2) - 2*omega_z_squared));
    double omega_minus = 1/2*(omega_0 - sqrt(pow(omega_0,2) - 2*omega_z_squared));

    double A_pluss = (v0 + omega_minus * x0)/(omega_minus - omega_pluss);
    double A_minus = -(v0 + omega_pluss * x0)/(omega_minus - omega_pluss);

    motion_r(0,0) = x0;
    motion_r(0,1) = 0.0;
    motion_r(0,2) = z0;

    for(int i = 1; i < steps; i++){
        t += dt;
        x = A_pluss*cos(-omega_pluss*t) + A_minus*cos(-omega_minus*t);
        y = A_pluss*sin(-omega_pluss*t) + A_minus*sin(-omega_minus*t);
        z = z0*cos(omega_z*t);
        motion_r.at(i,0) = x;
        motion_r.at(i,1) = y;
        motion_r.at(i,2) = z;
    }

    return motion_r;
}