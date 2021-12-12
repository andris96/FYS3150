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
 
  int M = 201;
  double h = 1.0/(M-1); 
  double dt = 2.5E-5;
  double T = 0.002;

  arma::vec x, y;

  x = arma::linspace(0,1,M);
  y = arma::linspace(0,1,M);
  
  double xc, yc, sx, sy, px, py;
  xc = 0.25;
  yc = 0.5;
  sx = 0.05;
  sy = 0.05;
  px = 200.0;
  py = 0.0;

  double v0 = 1.0E10;

  arma::cx_vec u = u_init(x,y,xc,yc,sx,sy,px,py);
  arma::mat V = V_config(2, v0, x, y);

  //simulate(u,V,dt,T,h);



  return 0;
}