#include "PenningTrap.hpp"
#include "Particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B0, double V0, double d) {
    B0_in = B0;
    V0_in = V0;
    d_in = d;

    k_e = 1.38935333*10e5; // Coulomb constant
    T = 9.64852558*10e1; // Magnetic field strength, Tesla
    V = 9.64852558*10e7; // Electric potential, Volt

    std::vector<Particle> particles;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in) {
    particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r) {
    // Get position values
    double x = r.at(0);
    double y = r.at(1);
    double z = r.at(2);

    // Compute the components of the gradient of V at $r
    double dx = -x*(V0_in/pow(d,2));
    double dy = -y*(V0_in/pow(d,2));
    double dz = 2*z*(V0_in/pow(d,2));

    // Compute the external electric field
    arma::vec E_field = arma::vec(3).fill(0.);
    E_field << -dx << -dy << -dz;
    return E_field;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r) {
    arma::vec B_field = arma::vec(3).fill(0.);
    B_field.at(2) = B0_in;
    return B_field;
}

// Force on particle_i from particle_j
// Assuming the force on p_i from p_j is govnered by the electric force only..
//  POSSIBLE ERRORS!!!!!
arma::vec PenningTrap::force_particle(int i, int j) {
    double q_i = particles.at(i).q();
    double q_j = particles.at(j).q();
    arma::vec r_i = particles.at(i).r();
    arma::vec r_j = particles.at(j).r();

    arma::vec force = arma::vec(3).fill(0.);
    force = k_e*q_i*q_j*(r_i - r_j)/(pow(arma::abs(r_i - r_j), 3));
    return force;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i) {
    // Initialize values
    arma::vec r_particle_i = particles.at(i).r();
    arma::vec v_particle_i = particles.at(i).v();
    arma::vec E_field = external_E_field(r_particle_i);
    arma::vec B_field = external_B_field(r_particle_i);
    double q_i = particles.at(i).q();

    // Compute the Lorentz force due to external fields
    arma::vec total_force = arma::vec(3).fill(0.);
    total_force = q_i*E_field + q_i*arma::cross(v_particle_i, B_field);
    return total_force;
}

// The total force on particle_i from the other particles
// POSSIBLE ERRORS!!!!
arma::vec PenningTrap::total_force_particles(int i) {
    //assert(particles.size > 1);
    arma::vec total_force = arma::vec(3).fill(0.);
    for (int j = 0; j < particles.size(); j++) {
        if (j != i) {
            total_force += force_particle(i, j); // k_e is baked into this method which gives severeal multiples of k_e..
        }
    }
    return total_force/pow(k_e, particles.size()-1); // fix this!!!! Had to remove multiples of k_e...
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i) {
    return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
// NOT COMPLETE!!!!
// POSSIBLE ERRORS!!!!
void PenningTrap::evolve_RK4(double dt){
    double k1, k2, k3, k4; 
    //Introducing a temporary particle to calculate the force when time and velocity has changed
    Particle pf; //pf = particle to calculate force
    for (int i = 0; i < particles.size(); i++){
        arma::vec v_i = particles.at(i).v();
        arma::vec r_i = particles.at(i).r();
        k1 = total_force(i)*dt; //excluding the mass, to lower number of calculations
        
        pf = particles.at(i); //pf = particle to calculate force
        pf.v() = v_i + k1/2; //changing v_i to v_i+(1/2)
        pf.r() = r_i + pf.v()*dt/2; //Changing the position with the new velocity

        k2 = total_force(pf)*dt; //Calculating k2 with the new force
        
        pf.v() = v_i + k2/2;
        pf.r() = r_i + pf.v()*dt/2;

        k3 = total_force(pf)*dt;

        pf.v() = v_i + k3;
        pf.r() = r_i + pf.v()*dt;

        k4 = total_force(pf)*dt;

        particles.at(i).v() = v_i + 1/6 * 1/particles.at(i).m() * (k1 + 2*k2 + 2*k3 + k4);
        particles.at(i).r() = particles.at(i).v()*dt;
    }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
    for (int i = 0; i < particles.size(); i++){
        particles.at(i).v() += total_force(i)*dt/particles.at(i).m();
        particles.at(i).r() += particles.at(i).v()*dt;
    }
}
