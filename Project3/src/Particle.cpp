#include "Particle.hpp"

// Constructors
Particle::Particle()
{
    q_ = 1;
    m_ = 1;
    r_ = arma::vec(3).fill(0.0);
    v_ = arma::vec(3).fill(0.0);
}

Particle::Particle(int q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
q_ = q_in;
m_ = m_in;
r_ = r_in;
v_ = v_in;
}

// Methods that returns the variables of the particle
// q = charge, m = mass, r = position, v = velocity
int Particle::q()
{
return q_;
}

double Particle::m()
{
return m_;
}

arma::vec Particle::r()
{
return r_;
}

arma::vec Particle::v()
{
return v_;
}

// Methods for changing the position and velocity of a particle
void Particle::set_r(arma::vec r_in) {
    r_ = r_in;
}

void Particle::set_v(arma::vec v_in) {
    r_ = v_in;
}

// Method for printing the position and velocity
void Particle::print_state() {
    std::cout << "Position: \n";
    r_.print();
    std::cout << "Velocity: \n";
    v_.print();
}