#include "functions.hpp"

/**
 * Takes indices i and j of a M x M matrix and convert them to a corresponding vector index k
*/
int convertk(int i, int j, int M){
    return (M*j + i);
}



/**
 * Constructs matrices A and B according to the Crank-Nicolson scheme in two dimensions
 * 
 * M is the number of points on x and y axis
 * s is the number of points on x and y axis excluding boundaries, so it is M-2
 * L is the number of rows and columns in A and B, which is the same as (M-2)*(M-2)
 * dt is the size of the timesteps
 * V is the potential matrix
 * A and B are the matrices we store the results in
 * a and b are the vectors containing the diagonal elements of A and B (not to be confused with the b = B*u in problem 3 and in main.cpp)
*/
void AB(int M, double h, double dt, arma::mat V, arma::sp_cx_mat &A, arma::sp_cx_mat &B){
    // Initiating
    int s = M-2;
    int L = s*s;
    int k = 0;
    std::complex<double> r = 1.i*dt/(2.*h*h);
    arma::cx_vec a(L, arma::fill::zeros);
    arma::cx_vec b(L, arma::fill::zeros);


    // Calculating the elements of a and b
    for(int i = 0; i < s; i++){
        for(int j = 0; j < s; j++){
            k = convertk(i,j,s);
            a(k) = 1. + 4.*r + 1.i * dt/2. * V(i,j);
            b(k) = 1. - 4.*r - 1.i * dt/2. * V(i,j);  
        }
    }

    // Setting the diagonal elements of A and B to a and b 
    for (int i = 0; i < L; i++){
        A(i,i) = a(i);
        B(i,i) = b(i);
    }

    // Setting up tridiagonal matrices
    A(0,0) = a(0);
    A(0,1) = -r;
    A(L-1,L-2) = -r;
    A(L-1,L-1) = a(L-1); 

    B(0,0) = b(0);
    B(0,1) = r;
    B(L-1,L-2) = r;
    B(L-1,L-1) = b(L-1); 

    for(int i=1; i<L-1; i++){
        A(i, i-1) = -r; 
        A(i, i+1) = -r;

        B(i, i-1) = r; 
        B(i, i+1) = r;
    }

    // Tridiagonal done, need to set every s element from the upper and lower diagonal to zero, in the case of M = 5, this would be every third element
    for (int i = s; i < L; i += s ){
        A(i,i-1) = 0.0;
        A(i-1,i) = 0.0;

        B(i,i-1) = 0.0;
        B(i-1,i) = 0.0;
    }

        
    // Setting up the rest r valued elements, these are the ones that goes diagonally directly next to the first submatrix with size (m-2)
    for(int i = 0; i < (L-s); i++){
        A(i,s+i) = -r;
        A(s+i,i) = -r;

        B(i,s+i) = r;
        B(s+i,i) = r;
    }

}

/**
 * Constructs the initial wave function as a gaussian wave packet
 * Takes x and y coordinates and calculates the wave function in that point with given parameters
 * xc and yc are the center coordinates of the packet 
 * sx and sy are the widths of the packets in x and y direction
 * px and py are the momentums in x and y direction
*/
arma::cx_vec u_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py){
    int M = x.size();
    int k = 0;
    arma::cx_vec u(M*M);
    for(int i = 0; i < M; i++){
        for(int j = 0; j < M; j++){
            k = convertk(i,j,M);
            if( (x(i) == 0.0) || (y(j) == 0.0) || (x(i) == 1.0) || (y(j) == 1.0) ){
                u(k) = 0.0;
            }
            else {
                u(k) = std::exp( -pow(x(i)-xc, 2)/(2.*sx*sx) - pow(y(j)-yc, 2)/(2.*sy*sy) + 1.i*px*(x(i)-xc) + 1.i*py*(y(j)-yc));
            }
        }
    }
    u = arma::normalise(u);
    return u;
}

/** 
 * In order to solve the matrix equations with A and B, we make an almost identical function where the boundaries are excluded 
 * Takes x and y coordinates and calculates the wave function in that point with given parameters
 * xc and yc are the center coordinates of the packet 
 * sx and sy are the widths of the packets in x and y direction
 * px and py are the momentums in x and y direction
 * 
 * This could be done in a different way by just altering the u_init to do both cases
*/
arma::cx_vec u_inner_init(arma::vec x, arma::vec y, double xc, double yc, double sx, double sy, double px, double py){
    int s = x.size()-2;
    int k = 0;
    arma::cx_vec u(s*s);
    for(int i = 0; i < s; i++){
        for(int j = 0; j < s; j++){
            k = convertk(i,j,s);
            u(k) = std::exp( -pow(x(i)-xc, 2)/(2.*sx*sx) - pow(y(j)-yc, 2)/(2.*sy*sy) + 1.i*px*(x(i)-xc) + 1.i*py*(y(j)-yc));
        }
    }
    u = arma::normalise(u);
    return u;
}

/**
 * Takes a M^2-dimentional vector and converts it to a M x M matrix by utilizing the function convertk
 */
arma::cx_mat vec_to_mat(arma::cx_vec u){
    int M = sqrt(u.size());
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

/**
 * Constructs the V matrix at x and y coordinates with either 1 2 or 3 slits
 * v0 is the potential of the barrier
 * The size of each slit is set to 0.05
 * And the thickness is set to 0.02
 */
arma::mat V_config(int slits, double v0, arma::vec x, arma::vec y){
    double slitsize = 0.05;
    double thickness = 0.02;
    arma::mat V = arma::mat(x.size(), y.size(), arma::fill::zeros);

    // Running a for-loop for each of the slit configurations 
    // could probably avoid code repetition by using more if-statements and having just a single for-loop
    // but this seems much more readable
    if (slits == 1){
        // construct a wall from y=0 to y=1 at x = 0.5
        // makes a slit with center at y = 0.5
        double slit_beginning = 0.5 - (slitsize/2.);
        double slit_end = 0.5 + (slitsize/2.);

        for(int i = 0; i < x.size(); i++){
            if ( x(i) > (0.5 - thickness) ){
                if ( x(i) < (0.5 + thickness)){
                    for(int j = 0; j < y.size(); j++){
                        V(i,j) = v0;
                        if(y(j) > slit_beginning){
                            if(y(j) < slit_end){
                                V(i,j) = 0.0;
                            }
                        }
                    }
                }     
            }
        }
    }

    if (slits == 2){
        // construct a wall from y=0 to y=1 at x = 0.5
        // makes 2 slits with center at y = 0.5+-slitsize
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

    if (slits == 3){
        // construct a wall from y=0 to y=1 at x = 0.5
        // makes 3 slits with center at y = 0.5, y = 0.5 - slitsize - slitsize/2, y =  0.5 + slitsize + slitsize/2
        double slit1_beginning = 0.5 - (slitsize/2.) - 2.*slitsize;
        double slit1_end = 0.5 - (slitsize/2.) - slitsize;
        double slit2_beginning = 0.5 - (slitsize/2.);
        double slit2_end = 0.5 + (slitsize/2.);
        double slit3_beginning = 0.5 + (slitsize/2.) + slitsize;
        double slit3_end = 0.5 + (slitsize/2.) + 2.*slitsize;

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
                        if (y(j) > slit3_beginning){
                            if(y(j) < slit3_end){
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
