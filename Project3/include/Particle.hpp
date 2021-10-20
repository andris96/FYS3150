#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <stdlib.h>
#include <armadillo>



class Particle {

    private:
    //Giving PenningTrap access to variables in Particle
    friend class PenningTrap;

    //Declearing variables for charge and mass
    int q_;
    double m_;
    arma::vec r_;
    arma::vec v_;

    public:
    //constructors
    Particle(); // Default, initiates with m=q=1 and r=v=0
    Particle(int c_in, double m_in, arma::vec r_in, arma::vec v_in);

    //Methods that returns the variables of the particle
    int q();
    double m();
    arma::vec r();
    arma::vec v();

    // Method for printing the position and velocity
    void print_state();

};

#endif