#include <iostream>
#include "Particle.hpp"
#include "PenningTrap.hpp"

int main()
{
    // Initating
    int q_ca = 1;
    double m_ca = 40.08; // [u]
    arma::vec r = arma::vec(3);
    arma::vec v = arma::vec(3);
    r << 1 << 1 << 1;
    v << 1 << 0 << 0;
    Particle p_ca = Particle(q_ca, m_ca, r, v);
    
    double tmax = 100; //100 micro seconds
    int steps = 10000;
    double dt = tmax/steps;

    double T = 9.64852558*10e1;
    double V = 9.64852558*10e7;

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 10*V; // Electric potential, Volt
    double d = 1000; // 1000 micrometers = 1 cm

    // Simulating a single particle in the trap
    //-----------------------------------------
    PenningTrap trap_euler = PenningTrap(B0, V0, d);
    PenningTrap trap_RK4 = PenningTrap(B0, V0, d);
    trap_euler.add_particle(p_ca);
    trap_RK4.add_particle(p_ca);

    arma::vec motion_z_euler(steps, arma::fill::zeros);
    arma::vec motion_z_RK4(steps, arma::fill::zeros);
    arma::vec time_interval = arma::linspace(0, tmax, steps);

    // Evolve the systems for $tmax micro seconds
    for(int i = 0; i < steps; i++){
        trap_euler.evolve_forward_Euler(dt);
        trap_RK4.evolve_RK4(dt);
        motion_z_euler.at(i) = trap_euler.get_particles().at(0).r().at(2);
        motion_z_RK4.at(i) = trap_RK4.get_particles().at(0).r().at(2);
    }

    motion_z_euler.save("motion_z_euler.txt", arma::raw_ascii);
    motion_z_RK4.save("motion_z_RK4.txt", arma::raw_ascii);
    time_interval.save("time_interval.txt", arma::raw_ascii);


    
    return 0;   
}