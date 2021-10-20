#include <iostream>
#include "/home/rajmann/FYS3150/Project3/include/Particle.hpp"

int main()
{
    // Intiate Particle object
    double q_in = 1;
    double m_in = 2.0;
    arma::vec v = arma::vec(3).fill(0.0);
    arma::vec r = arma::vec(3).fill(0.0);

    Particle p1 = Particle(q_in, m_in, v, r);

    std::cout << "Charge: " << p1.q() << '\n'
    << "Mass: " << p1.m() << '\n' << "Position: " << '\n'
    << p1.r() << '\n' << "Velocity: " << '\n' << p1.v() << std::endl;

    // Additional test
    Particle pp;

    return 0;
}