#include "Particle.hpp"

//Constructor
Particle::Particle(int q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
q_ = q_in;
m_ = m_in;
r_ = r_in;
v_ = v_in;
}

//Methods that returns the variables of the particle
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

