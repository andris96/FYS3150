#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <iomanip>
#include "functions.hpp"

using namespace std;



int main()
{
    int n = 100; //number of grid points
    //calculating values for x and u in problem 2
    vector<double> x = x_values(n);
    vector<double> exact_solution(n);
    
    //Creating a file which is called data.txt in order to save values of u and x
    ofstream myfile;  
    myfile.open ("data.txt"); 
    myfile << "x and u(x) \n";

  

    //time to calculate u as a function of x, and store values in the file data.txt
    for (int i = 0; i<n; i++ ){
    exact_solution[i] = u(x[i]); //calculating each u as function of x
    myfile << scientific << setprecision(3) << x[i] << " " << u(x[i]) << "\n";
    }

    myfile.close();  //closing the file
    
    

    return 0;
}
