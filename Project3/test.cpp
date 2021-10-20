#include <iostream>
#include "/home/rajmann/FYS3150/Project3/include/Particle.hpp"
#include "/home/rajmann/FYS3150/Project3/include/PenningTrap.hpp"

int main()
{
    arma::vec v = arma::vec(3).fill(0.0);
    arma::vec r = arma::vec(3).fill(0.0);
    Particle p1 = Particle(1, 2.0, v, r);
    
    double tmax = 100; //100 micro seconds
    int steps = 100; 
    double dt = tmax/steps;

    double T = 9.64852558*10e1; // Magnetic field strength, Tesla
    double V = 9.64852558*10e7; // Electric potential, Volt
    double d = 1000; //1000 micrometers = 1 cm

    PenningTrap trap_euler = PenningTrap(T,V,d);
    PenningTrap trap_RK4 = PenningTrap(T,V,d);
    trap_euler.add_particle(p1);
    trap_RK4.add_particle(p1);

    // Evolve the systems for $tmax micro seconds
    for(int i = 1; i < steps; i++){
        trap_euler.evolve_forward_Euler(dt);
        //trap_RK4.evolve_RK4(dt);
    }

    // Compare the state of the systems
    std::cout << "Trap evolved with Euler\n";
    trap_euler.print_states();

    std::cout << "Trap evolved with RK4\n";
    trap_RK4.print_states();

    return 0;   
}