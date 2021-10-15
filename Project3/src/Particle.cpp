#include "Particle.hpp"

//Constructor
Particle::Particle(int q_in, double m_in, arma::vec r_in, arma::vec v_in){
    q_ = q_in;
    m_ = m_in;
    r_ = r_in;
    v_ = v_in;

    //Methods that returns the variables of the particle
    int charge(){
        return q_;
    }

    double mass(){
        return m_;
    }

    arma::vec pos(){
        return r_;
    }

    arma::vec vel(){
        return v_;
    }

}