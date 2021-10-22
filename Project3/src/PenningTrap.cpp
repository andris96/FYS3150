#include "PenningTrap.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in) {
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
    particle_interaction_ = true;

    double k_e = 1.38935333*10e5; // Coulomb constant
    std::vector<Particle> particles;
}

// Getters
std::vector<Particle> PenningTrap::get_particles() {
    return particles;
}

bool PenningTrap::get_status_particle_interaction() {
    return particle_interaction_;
}

// Setters
void PenningTrap::enable_particle_interaction() {
    particle_interaction_ = true;
}

void PenningTrap::disable_particle_interaction() {
    particle_interaction_ = false;
}

// Print information about the state of the particles
void PenningTrap::print_states() {
    for (int i = 0; i < particles.size(); i++) {
        std::cout << "Particle nr. " << i << std::endl;
        particles.at(i).print_state();
        std::cout << "\n";
    }
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in) {
    particles.push_back(p_in);
}

// Check if a given position $r is within the boundaries of the trap
bool PenningTrap::is_within_trap(arma::vec r) {
    return (arma::norm(r) <= d_) ? true : false;
}


// Count the number of particles within the boundaries of the trap
int PenningTrap::count_particles() {
    int particle_count = 0;
    for (int i = 0; i < particles.size(); i++) {
        particle_count += (is_within_trap(particles.at(i).r())) ? 1 : 0;
    }
    return particle_count;
}


// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r) {
    if (is_within_trap(r)) {
        // Get position values
        double x = r.at(0);
        double y = r.at(1);
        double z = r.at(2);

        // Compute the components of the gradient of V at $r
        double dx = -x*(V0_/pow(d_,2));
        double dy = -y*(V0_/pow(d_,2));
        double dz = 2*z*(V0_/pow(d_,2));

        // Compute the external electric field
        arma::vec E_field = arma::vec(3).fill(0.);
        E_field << -dx << -dy << -dz;
        return E_field;
    }
    else {
        return arma::vec(3).fill(0.);
    }
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r) {
    if (is_within_trap(r)) {
        arma::vec B_field = arma::vec(3).fill(0.);
        B_field.at(2) = B0_;
        return B_field;
    }
    else {
        return arma::vec(3).fill(0.);
    }
}

// Force on particle_i from particle_j
// Assuming the force on p_i from p_j is govnered by the electric force only..
// POSSIBLE ERRORS!!!!!
arma::vec PenningTrap::force_particle(int i, int j) {
    if (particles.size() < 2) {
        return arma::vec(3).fill(0);
    }
    else{
    double q_i = particles.at(i).q();
    double q_j = particles.at(j).q();
    arma::vec r_i = particles.at(i).r();
    arma::vec r_j = particles.at(j).r();

    arma::vec force = arma::vec(3).fill(0.);
    force = k_e*q_i*q_j*(r_i - r_j)/(pow(arma::abs(r_i - r_j), 3));
    return force;
    }
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
arma::vec PenningTrap::total_force_particles(int i) {
    //No force if there is only 1 particle
    if (particles.size() < 2 || get_status_particle_interaction() == false) {
        return arma::vec(3).fill(0);
    }
    else{
    arma::vec total_force = arma::vec(3).fill(0.);

    for (int j = 0; j < particles.size(); j++) {
        if (j != i) {
            total_force += force_particle(i, j); 
        }
    }
    return total_force;
    }
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i) {
    return total_force_external(i) + total_force_particles(i);
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
// MAY NOT BE COMPLETE!!!!
// POSSIBLE ERRORS!!!!
void PenningTrap::evolve_RK4(double dt){
    arma::vec k1, k2, k3, k4; 
    for (int i = 0; i < particles.size(); i++){
        
        // Saving the state of particle i
        arma::vec v_i = particles.at(i).v_; 
        arma::vec r_i = particles.at(i).r_;

        k1 = total_force(i)*dt; // Excluding the mass, to lower number of calculations

        particles.at(i).v_ = v_i + k1/2; // Changing v_i to v_i+(1/2)
        particles.at(i).r_ = r_i + particles.at(i).v_*dt/2; // Changing the position with the new velocity
        k2 = total_force(i)*dt; // Calculating k2 with the new force
        
        particles.at(i).v_ = v_i + k2/2;
        particles.at(i).r_ = r_i + particles.at(i).v_*dt/2;
        k3 = total_force(i)*dt;

        particles.at(i).v_ = v_i + k3;
        particles.at(i).r_ = r_i + particles.at(i).v_*dt;
        k4 = total_force(i)*dt;

        // Taking mass into account here
        particles.at(i).v_ = v_i + 1./6. * (k1 + 2*k2 + 2*k3 + k4);
        particles.at(i).r_ = particles.at(i).v_*dt;

    }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
    for (int i = 0; i < particles.size(); i++){
        particles.at(i).v_ += total_force(i)*dt/particles.at(i).m_;
        particles.at(i).r_ += particles.at(i).v_*dt;
    }
}