#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "stdlib.h"
#include "math.h"
#include "functions.hpp"

using namespace std;



int main(){
    int n = 1000; //number of grid points
    //I managed to misunderstand the problem and tried to include the boundaries
    //which makes this variable actually m

    double h = 1/(double(n)-1); //stepsize

    vector<double> x = x_values(n);
    vector<double> func(n);

    for (int i = 0; i<n+1; i++){
        func[i] = f(x[i]);
    }

    //v has boundaries included, b is same size as v without boundaries
    //a and c are have 1 less element than b as they are the subdiagonal and superdiagonal
    //g is same size as b

    //the first element of a corresponds to the second element of b
    //which again corresponds to the third element of v and x

    /*
    a = [a_2,a_3,..,a_(n-1)]    size = n-3
    b = [b_1,b_2,..,b_(n-1)]    size = n-2
    c = [c_1,c_2,..,c_(n-2)]    size = n-3
    v = [v_0,v_1,..,v_n]
    x = [x_0,x_1,..,x_n]
    g = [g_1,g_2,..,g_(n-1)]    size = n-2
    */

    vector<double> a(n-3,-1);      
    vector<double> b(n-2,2);    
    vector<double> c(n-3,-1);    
    vector<double> g(n-2);        


    vector<double> b_tilde(n-2);
    b_tilde[0] = b[0]; 
    vector<double> g_tilde(n-2);
    g_tilde[0] = g[0];

    //v has the boundaries included and is therefore longer
    vector<double> v(n); 
    v[0] = 0;
    v[n] = 0;

    //calculating first and last element of g
    g[0] = pow(h,2)*func[0] + v[0];
    g[n-2] = pow(h,2)*func[n] + v[n];

    for (int i = 1; i<n-1; i++ ){
        g[i] = pow(h,2) * func[i];
    }
    
    

    //forward substitution
    
   
    for (int i = 1; i<n-1; i++ ){
        b_tilde[i] = b[i] - a[i-1]*c[i-1]/b_tilde[i-1]; 
        // the first iteration is then b_2-a_2*c_1/b~_1
        g_tilde[i] = g[i] - a[i-1]*g_tilde[i-1]/b_tilde[i-1];
    }


    //backward substitution

    v[n-1] = g_tilde[n-1]/b_tilde[n-1];
    for (int i = n-2; i>1; i-- ){
        v[i] = (g_tilde[i-1] - c[i-2]*v[i-1])/b_tilde[i-1];
    }
   
    
    ofstream myfile; 
    myfile.open ("v_x3.txt"); 
    myfile << "x and v(x) \n";

    for (int i = 0; i<n+1; i++ ){
    myfile << scientific << setprecision(5) << x[i] << " " << v[i] << "\n";
    }
    
    
    
    return 0;
}