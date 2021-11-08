#include "functions.hpp"


double degeneracy_energy(double epsilon, double total_states){
    return epsilon/total_states;
}

double degeneracy_magnetization(double m, double total_states){
    return m/total_states;
}

double probability_epsilon(double epsilon){}

double expectation_epsilon(double epsilon);

double probability_m(double m);

double expectation_m(double m);

double spesific_heat(double epsilon, double T);

double susceptibility(double m, double T);
