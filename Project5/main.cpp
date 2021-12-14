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
  /**
   * Defining values according to the task
   * M is the length of x and y vectors
   * A and B matrices are of length (L x L) 
  */
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
  double v0 = 1.0E10;

  xc = 0.25;
  yc = 0.5;
  sx = 0.05;
  sy = 0.2;
  px = 200.0;
  py = 0.0;
  

  arma::cx_vec u = u_init(x,y,xc,yc,sx,sy,px,py);
  arma::mat V = V_config(2, v0, x, y);
  arma::cx_mat U = vec_to_mat(u);
  arma::cx_vec u_inner = u_inner_init(x,y,xc,yc,sx,sy,px,py);
  arma::cx_mat U_inner = U.submat(1,1,M-2,M-2);
  arma::sp_cx_mat A(L,L);
  arma::sp_cx_mat B(L,L);
  arma::cx_vec b(L);
  arma::cx_cube U_t(M,M,tsteps);
  arma::cube p(M,M,tsteps);
  arma::vec p_sum(tsteps);
  

  AB(M,h,dt,V,A,B);


  for(int i = 0; i<tsteps; i++){
    b = B*u_inner;
    u_inner = arma::spsolve(A, b);
    U_inner = vec_to_mat(u_inner);
    U.submat(1,1,M-2,M-2) = U_inner;
    U_t.slice(i) = U;
    p.slice(i) = arma::real(arma::conj(U)*U);
    p_sum(i) = norm(U);
  }
 
  arma::cube U_real = arma::real(U_t);
  arma::cube U_imag = arma::real(U_t);

  p.save("p.bin", arma::raw_binary);
  p_sum.save("p_sum.bin", arma::raw_binary);
  U_t.save("U.bin", arma::raw_binary);
  U_real.save("U_real.bin", arma::raw_binary);
  U_imag.save("U_imag.bin", arma::raw_binary);

  // Problem 9
  // Finding the index where x = 0.8

  int i = 0;
  double index_x;
  while(x(i) < 0.8){
    index_x = i;
    i++;
  }

  arma::mat p_at_T = p.slice(tsteps-1);
  arma::vec p_given_x = p_at_T.col(index_x);

  arma::vec p_given_x_normalised = arma::normalise(p_given_x);

  p_given_x_normalised.save("p_given_x.bin", arma::raw_binary);




  return 0;
}