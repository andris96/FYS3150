#ifndef __TestIsingModel_hpp__
#define __TestIsingModel_hpp__

class TestIsingModel{
    private:
    //This class is a friend class of IsingModel
    //don't need anything that isn't already in IsingModel yet

    public:
    TestIsingModel::TestIsingModel(); // Empty constructor

    //calculating the different analytical solutions
    double TestIsingModel::Partition(double beta);
    double TestIsingModel::epsilon_squared_expectation(double beta);
    double TestIsingModel::m_squared_expectation(double beta);
    double TestIsingModel::Cv(double beta, double T);
    double TestIsingModel::khi(double beta, double T);

};

#endif