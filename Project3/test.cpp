#include <iostream>
#include "Particle.hpp"



int main()
{
    arma::vec v(3);
    arma::vec r(3);

    v.at(0) = 0.0;
    v.at(1) = 0.0;
    v.at(2) = 0.0;
    r.at(0) = 0.0;
    r.at(1) = 0.0;
    r.at(2) = 0.0;
    

    Particle p1 = Particle(1,2.0,v,r);

    std::cout << "Charge: " << p1.charge() << '\n'
    << "Mass: " << p1.mass() << '\n' << "Position: " << '\n'
    << p1.pos() << '\n' << "Velocity: " << '\n' << p1.vel() << std::endl;
}