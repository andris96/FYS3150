#include <iostream>
#include <string.h>

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
    Particle p_ca_1 = Particle(q_ca, m_ca, r, v);
    Particle p_ca_2 = Particle(q_ca, m_ca, -r, -v);

    
    double tmax = 100; //100 micro seconds
    int steps = 10000;
    double dt = tmax/steps;

    double T = 9.64852558*10e1;
    double V = 9.64852558*10e7;

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 10*V; // Electric potential, Volt
    double d = 1000; // 1000 micrometers = 1 cm

    // Simulating a single particle in the trap. Evolving the system
    // for $tmax micro seconds. Outputting z and t in .txt files for 
    // later plotting in python.
    //-----------------------------------------
    PenningTrap trap = PenningTrap(B0, V0, d);
    trap.add_particle(p_ca_1);

    arma::vec motion_z(steps, arma::fill::zeros);
    arma::vec time_interval = arma::linspace(0, tmax, steps);

    for(int i = 0; i < steps; i++){
        trap.evolve_RK4(dt);
        motion_z.at(i) = trap.get_particles().at(0).r().at(2);
    }

    motion_z.save("motion_z.txt", arma::raw_ascii);
    time_interval.save("time_interval.txt", arma::raw_ascii);

    // Simulating two particles in the trap. Evolving the system
    // for $tmax micro seconds.
    // Outputting x and y positions in .txt files for plotting 
    // in python. Simulating with and without particle interactions
    // enabled.
    //-----------------------------------------
    std::vector<std::string> interaction_mode(2);
    interaction_mode.at(0) = "with";
    interaction_mode.at(1) = "without";

    for (int mode = 0; mode <= 1; mode++) {

            PenningTrap trap = PenningTrap(B0, V0, d);
            if (interaction_mode.at(mode) == "without") {
                trap.disable_particle_interaction();
            }

            trap.add_particle(p_ca_1);
            arma::vec motion_x_1(steps, arma::fill::zeros);
            arma::vec motion_y_1(steps, arma::fill::zeros);

            trap.add_particle(p_ca_2);
            arma::vec motion_x_2(steps, arma::fill::zeros);
            arma::vec motion_y_2(steps, arma::fill::zeros);

            for(int i = 0; i < steps; i++){
                trap.evolve_RK4(dt);
                motion_x_1.at(i) = trap.get_particles().at(0).r().at(0);
                motion_x_2.at(i) = trap.get_particles().at(0).r().at(1);
                motion_y_1.at(i) = trap.get_particles().at(1).r().at(0);
                motion_y_2.at(i) = trap.get_particles().at(1).r().at(1);
            }
            motion_x_1.save("motion_x_1_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
            motion_y_1.save("motion_y_1_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
            motion_x_2.save("motion_x_2_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
            motion_y_2.save("motion_y_2_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
    }

    return 0;   
}