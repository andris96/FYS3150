#include <iostream>
#include <string.h>
#include <math.h>

#include "Particle.hpp"
#include "PenningTrap.hpp"

int main() {

    // Setting parameters for the simulation
    double T = 9.64852558*10e1; // Tesla, [u / e*(\mu*s)]
    double V = 9.64852558*10e7; // Volt, [u*(\mu*m)^2 / e*(\mu*s)^2]

    int q_ca = 1;
    double m_ca = 40.08; // [u]

    double B0 = 1*T;  // Magnetic field strength, Tesla
    double V0 = 0.0025*V; // Electric potential, Volt
    double d = 500; // 500 micrometers = 0.05 cm

    double tmax = 500; //500 micro seconds
    int t_steps = 100; 
    double dt = tmax/t_steps;

    double f = 0.1; // Amplitude <--- Set proper values!!!
    double omega_start = 1.0; // <--- Set proper values!!!
    double omega_end = 1.05;  // <--- Set proper values!!!
    double dOmega = 0.002;    // <--- Set proper values!!!
    int omega_steps = int((omega_end - omega_start) / dOmega);

    // Perform a fine-grained scan around in the frequency range [omega_start, omega_end]
    // ...
    // ...
    int numb_particles = 100;
    int numb_of_scans = 1;

    PenningTrap trap = PenningTrap(B0, V0, d);
    trap.enable_time_dependence();
    trap.set_f(f);

    std::vector<std::string> interaction_mode(2);
    interaction_mode.at(0) = "with";
    interaction_mode.at(1) = "without";

    // For each interaction mode
    for (int mode = 0; mode <= 1; mode++) {
        if (mode == 1) {
            trap.disable_particle_interaction();
        }

        arma::mat fractions(omega_steps, 2); 
        arma::vec result = arma::vec(2);

        // For each omega in the scan range
        for (int j = 0; j < omega_steps; j++) {
            double omega_V = omega_start + j*dOmega;

            // Taking the average of $numb_of_scans scans for each omega
            arma::vec fraction_omega(numb_of_scans);
            for (int s = 0; s < numb_of_scans; s++) {

                trap.fill_with_particles(q_ca, m_ca, numb_particles);
                trap.set_omega_V(omega_V);

                // Evolving the system for $tmax us
                for(int i = 0; i < t_steps; i++){
                    trap.evolve_RK4(dt);
                }

                fraction_omega.at(s) = trap.count_particles()/numb_particles;
                trap.reset();
            }

            double fraction = arma::mean(fraction_omega);
            result << omega_V << fraction;
            fractions.row(j)  = result.t();
        }
        fractions.save("fractions_fine_scan_" + interaction_mode.at(mode) + ".txt", arma::raw_ascii);
    }
    return 0;
}