#include <iostream>
#include "Particle.hpp"
#include "PenningTrap.hpp"



int main()
{
    arma::vec v(3);
    arma::vec r(3);

    v.at(0) = 0.0;
    v.at(1) = 0.0;
    v.at(2) = 0.0;
    r.at(0) = 0.0;
    r.at(1) = 0.0;
    r.at(2) = 0.0;

    Particle p1 = Particle(1,2.0,v,r);
    
    double tmax = 100; //100 microseconds
    int steps = 100; 
    double dt = tmax/steps;

    double T = 9.64852558*10e1; // Magnetic field strength, Tesla
    double V = 9.64852558*10e7; // Electric potential, Volt
    double d = 1000; //1000 micrometers = 1 cm

    PenningTrap trap = PenningTrap(T,V,d);
    trap.add_particle(p1);

    //not done
    for(int i = 0; i < steps; i++){
        trap.evolve_forward_Euler(dt);
        std::cout << p1.r().at(2) << std::endl;
    }

    
}