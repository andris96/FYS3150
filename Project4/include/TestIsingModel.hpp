#ifndef __TestIsingModel_hpp__
#define __TestIsingModel_hpp__

class TestIsingModel{
    private:
    //This class is a friend class of IsingModel
    //don't need anything that isn't already in IsingModel yet

    public:
    TestIsingModel(); // Empty constructor

    // Calculating the different analytical solutions
    double Partition(double beta);
    double epsilon_squared_expectation(double beta);
    double m_squared_expectation(double beta);
    double Cv(double beta, double T, int N);
    double khi(double beta, int N);

    // Methods for testing member methods of class IsingModel
    void generate_random_spin_config();

};

#endif