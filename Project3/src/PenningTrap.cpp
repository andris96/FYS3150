#include "PenningTrap.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in) {
    B0_ = B0_in;
    V0_ = V0_in;
    d_ = d_in;
    k_e_ = 1.38935333*10e5; // Coulomb constant, [u*(\mu*m)^2 / e*(\mu*s)^2]

    f_ = 1;
    omega_V_ = 1;

    particle_interaction_ = true;
    time_dependence_ = false;

    std::vector<Particle> particles_;
    total_time_ = 0; // Keep tracks of the over all time evolution of the system
}

// Reset the trap: clear all particles and reset the time
void PenningTrap::reset() {
    total_time_ = 0;
    particles_.clear();
}

// Getters
std::vector<Particle> PenningTrap::get_particles() {
    return particles_;
}
bool PenningTrap::get_status_particle_interaction() {
    return particle_interaction_;
}
bool PenningTrap::get_status_time_dependence() {
    return time_dependence_;
}
double PenningTrap::get_total_system_time() {
    return total_time_;
}

// Setters
void PenningTrap::enable_particle_interaction() {
    particle_interaction_ = true;
}
void PenningTrap::disable_particle_interaction() {
    particle_interaction_ = false;
}
void PenningTrap::enable_time_dependence() {
    time_dependence_ = true;
}
void PenningTrap::disable_time_dependence() {
    time_dependence_ = false;
}
void PenningTrap::set_f(double f_in) {
    f_ = f_in;
}
void PenningTrap::set_omega_V(double omega_V_in) {
    omega_V_ = omega_V_in;
}
// Print information about the state of the particles
void PenningTrap::print_states() {
    for (int i = 0; i < particles_.size(); i++) {
        std::cout << "Particle nr. " << i << std::endl;
        particles_.at(i).print_state();
        std::cout << "\n";
    }
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in) {
    particles_.push_back(p_in);
}

// Check if a given position $r is within the boundaries of the trap
bool PenningTrap::is_within_trap(arma::vec r) {
    return (arma::norm(r, 2) <= d_) ? true : false;
}

// Fill the trap with $n randomly initiated particles
void PenningTrap::fill_with_particles(int q, double m, int n) {
    for (int i = 0; i < n; i++) {
        arma::vec r = arma::vec(3).randn() * 0.1 * d_;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d_;  // random initial velocity
        Particle p(q, m, r, v);
        particles_.push_back(p);
    }
}

// Count the number of particles within the boundaries of the trap
int PenningTrap::count_particles() {
    int particle_count = 0;
    for (int i = 0; i < particles_.size(); i++) {
        particle_count += (is_within_trap(particles_.at(i).r())) ? 1 : 0;
    }
    return particle_count;
}

// External electric field at point r=(x,y,z)
// NOT COMPLETE!!
// Time dependence probably not correct..
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

        // Account for time-dependence
        return (get_status_time_dependence()) ? E_field*(1 + f_*cos(total_time_*omega_V_)) : E_field;//  <--- ?????
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
    if (particles_.size() < 2) {
        return arma::vec(3).fill(0);
    }
    else {
        double q_i = particles_.at(i).q();
        double q_j = particles_.at(j).q();
        arma::vec r_i = particles_.at(i).r();
        arma::vec r_j = particles_.at(j).r();

        arma::vec force = arma::vec(3).fill(0.);
        force = k_e_*q_i*q_j*(r_i - r_j)/(pow(arma::abs(r_i - r_j), 3));
        return force;
    }
}

// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i) {
    // Initialize values
    arma::vec r_particle_i = particles_.at(i).r();
    arma::vec v_particle_i = particles_.at(i).v();
    arma::vec E_field = external_E_field(r_particle_i);
    arma::vec B_field = external_B_field(r_particle_i);
    double q_i = particles_.at(i).q();

    // Compute the Lorentz force due to external fields
    arma::vec total_force = arma::vec(3).fill(0.);
    total_force = q_i*E_field + q_i*arma::cross(v_particle_i, B_field);
    return total_force;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i) {
    //No force if there is only 1 particle
    if (particles_.size() < 2 || get_status_particle_interaction() == false) {
        return arma::vec(3).fill(0);
    }
    else{
    arma::vec total_force = arma::vec(3).fill(0.);

    for (int j = 0; j < particles_.size(); j++) {
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

// Evolve the system one time step (dt) using Runge-Kutta 4th order.
// The time-depedence of the electric field is "baked into" the
// method external_E_field() and kept track of by the class, hence
// RK4 only needs to update the total time elapsed.
//
// MAY NOT BE COMPLETE!!!!
// POSSIBLE ERRORS!!!!
void PenningTrap::evolve_RK4(double dt){
    arma::vec k1, k2, k3, k4; 
    for (int i = 0; i < particles_.size(); i++){
        // Mass of particle
        double m = particles_.at(i).m_;
        // Saving the state of particle i
        arma::vec v_i = particles_.at(i).v_; 
        arma::vec r_i = particles_.at(i).r_;

        k1 = total_force(i)/m*dt; // Excluding the mass, to lower number of calculations

        particles_.at(i).v_ = v_i + k1/2; // Changing v_i to v_i+(1/2)
        particles_.at(i).r_ = r_i + particles_.at(i).v_*dt/2; // Changing the position with the new velocity
        k2 = total_force(i)/m*dt; // Calculating k2 with the new force
        
        particles_.at(i).v_ = v_i + k2/2;
        particles_.at(i).r_ = r_i + particles_.at(i).v_*dt/2;
        k3 = total_force(i)/m*dt;

        particles_.at(i).v_ = v_i + k3;
        particles_.at(i).r_ = r_i + particles_.at(i).v_*dt;
        k4 = total_force(i)/m*dt;

        // Taking mass into account here
        particles_.at(i).v_ = v_i + 1./6. * (k1 + 2*k2 + 2*k3 + k4);
        particles_.at(i).r_ = particles_.at(i).v_*dt;
    }

    // Keeping track of the time-evolution
    total_time_ += dt;
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt){
    for (int i = 0; i < particles_.size(); i++){
        particles_.at(i).v_ += total_force(i)*dt/particles_.at(i).m_;
        particles_.at(i).r_ += particles_.at(i).v_*dt;
    }
    total_time_ += dt;
}