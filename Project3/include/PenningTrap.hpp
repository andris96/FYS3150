#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include "Particle.hpp"
#include <stdlib.h>

// redundant..
#include <vector>
#include <armadillo>

class PenningTrap {

    protected:
    double B0_, V0_, d_;
    double k_e_;
    double f_, omega_V_;
    bool particle_interaction_;
    bool time_dependence_;
    std::vector<Particle> particles_;
    double total_time_;
    
    public:
    // Constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Reset the trap: clear all particles and reset the time
    void reset();

    // Getters
    std::vector<Particle> get_particles();
    bool get_status_particle_interaction();
    bool get_status_time_dependence();
    double get_total_system_time();
    
    // Setters
    void enable_particle_interaction();
    void disable_particle_interaction();
    void enable_time_dependence();
    void disable_time_dependence();
    void set_f(double f_in);
    void set_omega_V(double omega_V_in);

    // Print information about the state of the particles
    void print_states();

    // Add a particle to the trap
    void add_particle(Particle p_in);

    // Fill the trap with $n randomly initiated particles
    void fill_with_particles(int q, double m, int n);

    // Check if a given position $r is within the boundaries of the trap
    bool is_within_trap(arma::vec r);

    // Count the number of particles within the boundaries of the trap
    int count_particles();

    // External electric field at point r=(x,y,z)
    arma::vec external_E_field(arma::vec r);  

    // External magnetic field at point r=(x,y,z)
    arma::vec external_B_field(arma::vec r);  

    // Force on particle_i from particle_j
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the external fields
    arma::vec total_force_external(int i);

    // The total force on particle_i from the other particles
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles
    arma::vec total_force(int i);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order
    void evolve_RK4(double dt);

    // Evolve the system one time step (dt) using Forward Euler
    void evolve_forward_Euler(double dt);
    
};

#endif