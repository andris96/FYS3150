#include <iostream>
#include "Particle.hpp"
#include "PenningTrap.hpp"

int main()
{
    arma::vec v = arma::vec(3).fill(1);
    arma::vec r = arma::vec(3).fill(1);
    double m_ca = 40.08; // [u]
    int q_ca = 1;
    Particle p1 = Particle(q_ca, m_ca, v, r);
    
    double tmax = 100; //100 micro seconds
    int steps = 100; 
    double dt = tmax/steps;

    double T = 9.64852558*10e1;
    double V = 9.64852558*10e7;

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 10*V; // Electric potential, Volt
    double d = 1000; // 1000 micrometers = 1 cm

    PenningTrap trap_euler = PenningTrap(B0, V0, d);
    PenningTrap trap_RK4 = PenningTrap(B0, V0, d);
    trap_euler.add_particle(p1);
    trap_RK4.add_particle(p1);

    arma::vec motion_z(steps, arma::fill::zeros);
    arma::vec time_interval = arma::linspace(0,tmax,steps);

    // Evolve the systems for $tmax micro seconds
    for(int i = 0; i < steps; i++){
        trap_euler.evolve_forward_Euler(dt);
        motion_z(i) = p1.r().at(2);
    }

    motion_z.save("motion_z.txt", arma::raw_ascii);
    time_interval.save("time_interval.txt", arma::raw_ascii);
   
    return 0;   
}