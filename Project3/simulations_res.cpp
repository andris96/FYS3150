#include <iostream>
#include <string.h>
#include <math.h>

#include "Particle.hpp"
#include "PenningTrap.hpp"

#include <chrono> // tmp

// To build: 
// g++ simulations_res.cpp src/Particle.cpp src/PenningTrap.cpp src/utils.cpp
// -I include -o simulations_res.exe -larmadillo
//
// To run: ./simulations_res.exe


int main() {

    auto start = std::chrono::high_resolution_clock::now(); // tmp

    // Setting parameters for the simulation
    //double k_e = 1.38935333*10e5; // Coulomb constant, [u*(\mu*m)^2 / e*(\mu*s)^2]
    double T = 9.64852558*10e1; // Tesla, [u / e*(\mu*s)]
    double V = 9.64852558*10e7; // Volt, [u*(\mu*m)^2 / e*(\mu*s)^2]

    int q_ca = 1;
    double m_ca = 40.08; // [u]

    // Tmp notes: 
    // tmax=50, t_steps=1000 -> ca. 4 min 24 sec     fractions = 1.0 for all f
    // tmax=50, t_steps=500  -> ca. 2 min 14 sec     fractions = 0.0 for all f
    // tmax=500, t_steps=10000  -> ca. 31 min           fractions 1

    double tmax = 500; //500 micro seconds
    int t_steps = 10000; 
    double dt = tmax/t_steps;

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 0.0025*V; // Electric potential, Volt
    double d = 500; // 500 micrometers = 0.05 cm

    // Instantiating trap with particle interactions disabled
    PenningTrap trap = PenningTrap(B0, V0, d);
    trap.disable_particle_interaction();

    // Initiating with 100 randomly initialized Ca-particles. For each of
    // the amplitudes f = {0.1, 0.4, 0.7}, the fraction of particles still 
    // inside the trap after 500 microseconds is recorded as a function 
    // of the angular frequency \omgega_V \in [0.2, 2.5] MHz. 
    // The results are exported to .txt files for later plotting in 
    // python.
    //--------------------------------------------------------------
    std::vector<double> f_vec{0.1, 0.4, 0.7};
    
    double omega_start = 0.2;
    double omega_end = 2.5;
    double dOmega = 0.02;
    int omega_steps = int((omega_end - omega_start) / dOmega);

    int numb_particles = 100; // <---

    for (int k = 0; k < f_vec.size(); k++) {
        // For storing results, row-format: [omega_V fraction]
        arma::mat fractions(omega_steps, 2); 
        arma::vec result = arma::vec(2);

        double f = f_vec.at(k);
        for (int j = 0; j < omega_steps; j++) {
            double omega_V = omega_start + j*dOmega;

            trap.fill_with_particles(q_ca, m_ca, numb_particles);
            trap.set_f(f);
            trap.set_omega_V(omega_V);

            // Evolving the system for 500 us
            for(int i = 0; i < t_steps; i++){
                trap.evolve_RK4(dt);
            }

            double fraction = trap.count_particles()/numb_particles;
            result << omega_V << fraction;
            fractions.row(j)  = result.t();

            trap.reset();
        }

        fractions.save("fractions_f_" + std::to_string(f).substr(0, 4) + ".txt", arma::raw_ascii);
    }


    auto stop = std::chrono::high_resolution_clock::now(); // tmp
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}