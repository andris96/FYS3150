  #include "functions.hpp" 



/**
 * 
 * To build:
 * g++ main.cpp src/functions.cpp -I include -o main.exe -larmadillo
 * 
 * To run:
 * ./main.exe
*/



int main()
{

  int M = 5;
  int L = (M-2)*(M-2);

  
  double dt = 0.01;
  double h = 0.01;
  arma::mat V(L, L, arma::fill::randu);
  arma::cx_mat A(L,L, arma::fill::zeros);
  arma::cx_mat B(L,L, arma::fill::zeros);

  AB(M,h,dt,V,A,B);
  
  std::cout << std::setw(3) << A << std::endl;
  std::cout << std::setw(3) << B << std::endl;
  


    return 0;
}