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

  int M = 100;
  int L = (M-2)*(M-2);
  double t = 1;
  double t0 = 0;
  int tsteps = 100;
  int xysteps = 100;
  
  double dt = (t-t0)/tsteps;

  double h = 0.01;
  arma::mat V(L, L, arma::fill::zeros);
  arma::sp_cx_mat A(L,L);
  arma::sp_cx_mat B(L,L);

  AB(M,h,dt,V,A,B);

  arma::vec xy_values = arma::linspace(0,1,xysteps);
  arma::cx_vec x(L);
  arma::cx_vec y(L);
  for(int i; i < L; i++){
    x(i) = xy_values(i);
    y(i) = xy_values(i);
  }

  arma::cx_mat u(L,L);
  arma::cx_vec b(L);

  double xc, yc, sx, sy, px, py;
  xc = 0.;
  yc = 0.;
  sx = 0.;
  sy = 0.;
  px = 0.;
  py = 0.;
  

  

  
  


  //std::cout << std::setw(3) << "real A \n" << arma::real(A) << std::endl;
  //std::cout << std::setw(3) << "imag A \n" << arma::imag(A) << std::endl;
  //std::cout << std::setw(3) << "real B \n" << arma::real(B) << std::endl;
  //std::cout << std::setw(3) << "imag B \n" << arma::imag(B) << std::endl;

  


    return 0;
}