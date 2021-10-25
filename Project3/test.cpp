#include <iostream>

#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "TestPenningTrap.hpp"

// To compile:
// g++ test.cpp src/PenningTrap.cpp src/Particle.cpp -I include -o test.exe -larmadillo


int main()
{
    // int q_ca = 1;
    // double m_ca = 40.08; // [u]
    // arma::vec v = arma::vec(3);
    // arma::vec r = arma::vec(3);
    // v << 1 << 0 << 0;
    // r << 1 << 1 << 1;
    // Particle p1 = Particle(q_ca, m_ca, r, v);
    
    // double tmax = 100; //100 micro seconds
    // int steps = 10000; 
    // double dt = tmax/steps;

    // double T = 9.64852558*10e1;
    // double V = 9.64852558*10e7;

    // double B0 = 1*T;  // Magnetic field strength, Tesla
    // double V0 = 10*V; // Electric potential, Volt
    // double d = 10e4; // 10000 micrometers = 1 cm

    // PenningTrap trap_euler = PenningTrap(B0, V0, d);
    // PenningTrap trap_RK4 = PenningTrap(B0, V0, d);
    // trap_euler.add_particle(p1);
    // trap_RK4.add_particle(p1);

    // // Evolve the systems for $tmax micro seconds
    // for(int i = 0; i < steps; i++){
    //     trap_euler.evolve_forward_Euler(dt);
    //     trap_RK4.evolve_RK4(dt);
    // }

    // // Compare the states of the particles in both systems
    // std::cout << "Time evolution: " << tmax << " us " << "with " << steps << " steps\n\n"
    // << "Trap evolved with Euler\n";
    // trap_euler.print_states();

    // std::cout << "Trap evolved with RK4 at\n";
    // trap_RK4.print_states();

    // std::cout << "External E.field:" << '\n' << trap_euler.external_E_field(p1.r());
    // std::cout << "External B.field:" << '\n' << trap_euler.external_B_field(p1.r());
    // std::cout << "Total external forces:" << '\n' << trap_euler.total_force_external(0);
    // std::cout << "Total force from other particles:" << '\n' << trap_euler.total_force_particles(0);
    // std::cout << "Total force:" <<  '\n' << trap_euler.total_force(0);
    // std::cout << "Position of particle after iterations:" <<  '\n' << p1.r();

    // double d = 1;
    // arma::vec r = arma::vec(3).fill(0.9*d);

    // assert(arma::norm(r, 1) > 1);

    // Testing PenningTrap class
    TestPenningTrap test_trap = TestPenningTrap();
    test_trap.runAllTests();

    return 0;
}