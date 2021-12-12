#include "functions.hpp"

/**
 * Takes indices i and j of a M^2 x M^2 matrix and convert them to a corresponding vector index k
*/
int convertk(int i, int j, int M){
    return (M*j + i);
}



/**
 * Constructs matrices A and B according to the Crank-Nicolson scheme in two dimensions
 * 
*/
void AB(int M, double h, double dt, arma::mat V, arma::sp_cx_mat &A, arma::sp_cx_mat &B){
    
    int L = pow((M-2),2);
    std::complex<double> r = 1.i*dt/(2.*h*h);

    arma::cx_vec a(L, arma::fill::zeros);
    arma::cx_vec b(L, arma::fill::zeros);

    std::cout << "V = " << V.size() << std::endl;
    std::cout << L << std::endl;
    std::cout << A.size() << std::endl;

    int k = 0;

    //Calculating the elements of a and b
    for(int i = 0; i < M-2; i++){
        for(int j = 0; j < M-2; j++){
            if (i == j){
                k = convertk(i,j,M-2);
                std::cout << k << std::endl;
                a(k) = 1. + 4.*r + 1.i * dt/2. * V(i,j);
                b(k) = 1. - 4.*r - 1.i * dt/2. * V(i,j);
                A(i,j) = a(k);
                B(i,j) = b(k);
            }
        }
    }

    A(0,0) = a(0);
    A(0,1) = -r;
    A(L-1,L-2) = -r;
    A(L-1,L-1) = a(L-1); 

    B(0,0) = b(0);
    B(0,1) = r;
    B(L-1,L-2) = r;
    B(L-1,L-1) = b(L-1); 


    // Setting up a tridiagonal matrix
    for(int i=1; i<L-1; i++){
        A(i, i-1) = -r; 
        A(i, i+1) = -r;

        B(i, i-1) = r; 
        B(i, i+1) = r;
    }

        
    // Setting up the rest r valued elements, these are the ones that goes diagonally directly next to the first submatrix with size (m-2)
    int s = sqrt(L);
    for(int i = 0; i < (L-s); i++){
        A(i,s+i) = -r;
        A(s+i,i) = -r;

        B(i,s+i) = r;
        B(s+i,i) = r;
    }

}


arma::cx_vec u_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py){
    int M = x.size();
    int k = 0;
    arma::cx_vec u(M*M);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            k = convertk(i,j,M);
            if( (x(i) == 0.0) || (y(j) == 0.0) ){
                u(k) = 0.0;
            }
            else {
                u(k) = std::exp( -pow(x(i)-xc, 2)/(2.*sx*sx) - pow(y(i)-yc, 2)/(2.*sy*sy) + 1.i*px*(x(i)-xc) + 1.i*py*(y(i)-yc));
            }
        }
    }
    u = arma::normalise(u);
    return u;
}

arma::cx_mat vec_to_mat(arma::cx_vec u){
    int M = u.size();
    arma::cx_mat U(M,M);
    int k = 0;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            k = convertk(i,j,M);
            U(i,j) = u(k); 
        }
    }
    return U;
}


arma::mat V_config(int slits, double v0, arma::vec x, arma::vec y){
    double slitsize = 0.05;
    double thickness = 0.02;
    arma::mat V = arma::mat(x.size(), y.size(), arma::fill::zeros);
    if (slits == 2){
        // construct a wall from y=0 to y=1 at x = 0.5
        // make slits at 0.5+-0.025 that are 0.02 thick
        double slit1_beginning = 0.5 - (slitsize/2.) - slitsize;
        double slit1_end = 0.5 - (slitsize/2.);
        double slit2_beginning = 0.5 + (slitsize/2.);
        double slit2_end = 0.5 + (slitsize/2.) + slitsize;

        for(int i = 0; i < x.size(); i++){
            if ( x(i) > (0.5 - thickness) ){
                if ( x(i) < (0.5 + thickness)){

                    for(int j = 0; j < y.size(); j++){
                        V(i,j) = v0;
                        if(y(j) > slit1_beginning){
                            if(y(j) < slit1_end){
                                V(i,j) = 0.0;
                            }
                        }
                        if (y(j) > slit2_beginning){
                            if(y(j) < slit2_end){
                                V(i,j) = 0.0;
                            }
                        }
                    }
                }     
            }
        }
    }

    return V;
}

void simulate(arma::cx_vec u, arma::mat V, double dt, double T, double h){
    
    int tsteps = T/dt;
    int M = u.size();
    int L = (M-2)*(M-2);
    arma::cx_vec b(L);
    arma::sp_cx_mat A(L,L);
    arma::sp_cx_mat B(L,L);
    arma::cx_cube U_t(L,L,tsteps);

    AB(M,h,dt,V,A,B);
    
    for(int i = 0; i<tsteps; i++){
        b = B*u;
        u = arma::spsolve(A, b, "lapack");
        U_t.slice(i) = vec_to_mat(u);
    }
    
    U_t.save("U.bin");
}