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
  int L = (M-2)*(M-2);
  double h = 1.0/(M-1); 
  double dt = 2.5E-5;
  double T = 0.002;
  int tsteps = T/dt;

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
  arma::cx_mat U = vec_to_mat(u);
  arma::sp_cx_mat A(L,L);
  arma::sp_cx_mat B(L,L);
  arma::cx_cube U_t(M,M,tsteps);
  
  //std::cout << u.size() << std::endl;

  AB(M,h,dt,V,A,B);


  arma::cx_vec u_inner = u_inner_init(x,y,xc,yc,sx,sy,px,py);
  arma::cx_mat U_inner = U.submat(1,1,M-2,M-2);
  arma::cx_vec b(L);
  int k = 0;
  

  //std::cout << arma::real(U_inner) << std::endl;

  for(int i = 0; i<tsteps; i++){
    b = B*u_inner;
    u_inner = arma::spsolve(A, b);
    U_inner = vec_to_mat(u_inner);
    U.submat(1,1,M-2,M-2) = U_inner;
    U_t.slice(i) = U;
  }

  U_t.save("U.bin");


  return 0;
}