#include <iostream>
#include <string.h>
#include <math.h>

#include "Particle.hpp"
#include "PenningTrap.hpp"
#include "utils.hpp"

// To compile: 
// g++ simulations.cpp src/Particle.cpp src/PenningTrap.cpp src/utils.cpp
// -I include -o simulations.exe -larmadillo
// To run: ./simulations.exe

int main()
{
    // Initiating
    //double k_e = 1.38935333*10e5; // Coulomb constant, [u*(\mu*m)^2 / e*(\mu*s)^2]
    double T = 9.64852558*10e1; // Tesla, [u / e*(\mu*s)]
    double V = 9.64852558*10e7; // Volt, [u*(\mu*m)^2 / e*(\mu*s)^2]

    int q_ca = 1;
    double m_ca = 40.08; // [u]
    arma::vec r = arma::vec(3);
    arma::vec v = arma::vec(3);
    r << -.1 << 0.0 << .1;
    v << 0.0 << -.1 << 0.0;

    //parameters for the analytical solution
    double x0 = r(0);
    double z0 = r(2);
    double v0 = v(1);

    Particle p_ca_1 = Particle(q_ca, m_ca, r, v);
    Particle p_ca_2 = Particle(q_ca, m_ca, -r, -v);

    double tmax = 100; //100 micro seconds
    int steps = 1000;
    double dt = tmax/steps;

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 10*V; // Electric potential, Volt
    double d = 1000; // 1000 micrometers = 1 cm

    // Simulating a single particle in the trap. Evolving the system
    // for $tmax micro seconds. Outputting z and t in .txt files for 
    // later plotting in python.
    //--------------------------------------------------------------
    PenningTrap trap = PenningTrap(B0, V0, d);
    trap.add_particle(p_ca_1);

    arma::vec motion_z(steps, arma::fill::zeros);
    arma::mat motion_analytical = solve_analytical_1p(v0, x0, z0, tmax, steps, q_ca, B0, V0, m_ca, d);
    arma::vec motion_z_analytical(steps, arma::fill::zeros);
    arma::vec time_interval = arma::linspace(0, tmax, steps);

    for(int i = 0; i < steps; i++){
        trap.evolve_RK4(dt);
        motion_z.at(i) = trap.get_particles().at(0).r().at(2);
        motion_z_analytical(i) = motion_analytical(i,2);
    }

    


    motion_z.save("motion_z_RK4.txt", arma::raw_ascii);
    motion_z_analytical.save("motion_z_analytical.txt", arma::raw_ascii);
    time_interval.save("time_interval.txt", arma::raw_ascii);

    // Simulating two particles in the trap. Evolving the system
    // for $tmax micro seconds.
    // Outputting positions and velocities in .txt files for plotting 
    // in python. Simulating with and without particle interactions
    // enabled. 
    //(Generates what's needed for tasks 2-4 in problem 9).
    //--------------------------------------------------------------
    std::vector<std::string> interaction_mode(2);
    interaction_mode.at(0) = "with";
    interaction_mode.at(1) = "without";

    for (int mode = 0; mode <= 1; mode++) {
        PenningTrap trap = PenningTrap(B0, V0, d);
        if (interaction_mode.at(mode) == "without") {
            trap.disable_particle_interaction();
        }

        trap.add_particle(p_ca_1);
        arma::mat motion_r_1(steps, 3, arma::fill::zeros);
        arma::mat motion_v_1(steps, 3, arma::fill::zeros);

        trap.add_particle(p_ca_2);
        arma::mat motion_r_2(steps, 3, arma::fill::zeros);
        arma::mat motion_v_2(steps, 3, arma::fill::zeros);

        for(int i = 0; i < steps; i++){
            trap.evolve_RK4(dt);
            motion_r_1.row(i) = trap.get_particles().at(0).r().t(); // .t() : transposing col -> row
            motion_v_1.row(i) = trap.get_particles().at(0).v().t();
            motion_r_2.row(i) = trap.get_particles().at(1).r().t();
            motion_v_2.row(i) = trap.get_particles().at(1).v().t();
        }

        motion_r_1.save("motion_r_1_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
        motion_v_1.save("motion_v_1_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
        motion_r_2.save("motion_r_2_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
        motion_v_2.save("motion_v_2_" + interaction_mode.at(mode) + "_interactions.txt", arma::raw_ascii);
    }

    // Once again simulating a single particle in the trap. 
    // Evolving the system for 50 micro seconds in five different step 
    // lenghts in order to compare the numerical with analytical.
    // In addtion, the convergence rates for Euelr and RK4 is
    // computed.
    //(Generates what's needed for tasks 5-6 in problem 9).
    //--------------------------------------------------------------
    tmax = 50;
    arma::vec steps_vec = arma::vec(5);
    steps_vec << 50 << 100 << 500 << 1000 << 5000;

    std::vector<std::string> methods(2);
    methods.at(0) = "rk4";
    methods.at(1) = "euler";

    std::vector<double> delta_max_euler;
    std::vector<double> delta_max_rk4;
    std::vector<double> h_vec;

    for (int method = 0; method <= 1; method++) {
        for (int s = 0; s < steps_vec.size(); s++) {
            PenningTrap trap = PenningTrap(B0, V0, d);
            trap.add_particle(p_ca_1);

            arma::mat motion_r = arma::mat(steps_vec.at(s), 3, arma::fill::zeros);
            arma::mat motion_r_analytical = solve_analytical_1p(v0, x0, z0, tmax, 
                                                                steps_vec.at(s), q_ca, B0, V0, m_ca, d);
            arma::mat motion_r_diff = arma::mat(steps_vec.at(s), 3, arma::fill::zeros);

            // Evolving the system
            dt = tmax/steps_vec.at(s);
            for(int i = 0; i < steps_vec.at(s); i++){
                if (methods.at(method) == "rk4") {
                    trap.evolve_RK4(dt);
                }
                else {
                    trap.evolve_forward_Euler(dt);
                }
                motion_r.row(i) = trap.get_particles().at(0).r().t();
            }

            // Storing results for computing convergence rates later
            motion_r_diff = arma::abs(motion_r_analytical - motion_r);
            double delta_max = motion_r_diff.max();
            if (methods.at(method) == "rk4") {
                delta_max_rk4.push_back(delta_max);
                h_vec.push_back(dt);
            }
            else {
                delta_max_euler.push_back(delta_max);
            }

            // Saving for later plotting in python
            motion_r.save("motion_r_h" + std::to_string(s) + "_" + methods.at(method) + ".txt", arma::raw_ascii);
            motion_r_analytical.save("motion_r_h" + std::to_string(s) + "_analytical.txt", arma::raw_ascii);
        }
    }

    // Computing convergence rates and printing to terminal
    double convergence_rate_euler;
    double convergence_rate_rk4;
    for (int k = 1; k < h_vec.size(); k++) {
        convergence_rate_euler += log(delta_max_euler.at(k)/delta_max_euler.at(k-1)) / log(h_vec.at(k)/h_vec.at(k-1));
        convergence_rate_rk4 += log(delta_max_rk4.at(k)/delta_max_rk4.at(k-1)) / log(h_vec.at(k)/h_vec.at(k-1));
    }
    convergence_rate_euler /= 4.;
    convergence_rate_rk4 /= 4.;

    std::cout << "Convergence rate, Euler: " << convergence_rate_euler
    << "\nConvergence rate, RK4: " << convergence_rate_rk4 << std::endl;

    return 0;   
}