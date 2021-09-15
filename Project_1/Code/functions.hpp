#include <vector>
#include "math.h"

using namespace std;

//I want values in the range 0 to 1, where n is the amount of steps
//Making a function that returns a vector which contains all values of x
vector<double> x_values(int n){ 

    vector<double> x(n);           //declaring vector of size n
    double stepsize = 1/(double(n) -1);    //calculating stepsize, xmax-xmin/steps = 1/steps

    for (int i = 0; i<n; i++ ){
        x[i] += i*stepsize;      //assigning x-values
    }
    
    return x;

}

double f(double x){
    return 100*exp(-10*x);
}

double u(double x){
    double solution;
    if(x == 1 || x == 0){
        solution = 0;   //setting the boundary u(0) = 0, u(1) = 0
    }
    else{
    solution = 1-(1 - exp(-10))*x - exp(-10*x); //The function defined in problem 1
    }
    return solution;
}