#include <armadillo>
#include "IsingModel.hpp"
#include "TestIsingModel.hpp"

//Empty constructor
TestIsingModel::TestIsingModel() = default;

//Don't know how to get the variables from TestIsingModel
//Could save computation by calculating exp(8*beta) and saving that value

double TestIsingModel::Partition(double beta){
    double Z = 4*exp(8*beta) + 12;
    return Z;
}

double TestIsingModel::epsilon_squared_expectation(double beta){
    double eps_sq = 16*exp(8*beta)/Partition(beta); 
    //could make this less computational by having Partition(beta) saved as a constant somehow
    return eps_sq;
}

double TestIsingModel::m_squared_expectation(double beta){
    double m_sq = (2*exp(8*beta)+4)/Partition(beta);
    return m_sq;
}

double TestIsingModel::Cv(double beta, double T){
    double kB = 1; // should get this from IsingModel somehow
    double Cv = 1/(kB*T*T) * 16*exp(8*beta);
    return Cv;
}

double TestIsingModel::khi(double beta, double T){
    double khi = beta*(2*(exp(8*beta) + 1)/Partition(beta));
    return khi;
}