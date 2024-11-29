#include "Include/psi-boil.h"
#include <fstream>
#include <vector>

//Domian
const real L_l = 0.015; //length of the fluid domain
const real L_s = 0.0025; //length of the solid domain
double x_in = 0.015;    //location of the interface
const real T = 10;     //total run time of the simulation

const int nx_l = 20;   //grid points on fluid domain
const int nx_s = 20;   //grid points on solid domain
const int nt = 100;    //iterations

const double pi = acos(-1.0);

//Properties of fluid and solid

const double lam_s = 401.0; // thermal conductivity of copper
const double lam_f = 0.56; // thermal conductivity of Water
const double Cp_f = 4184; // water J/kg·K
const double Cp_s = 390; //  steel J/kg·K
const double rho_l = 1000;
const double rho_s = 8960;


const double T_bc = 373.0; // 100*C
const double q = -1000;//

/******************************************************************************/

 double **createZeroMatrix(int n, int m)
 {
    double **result = (double **) malloc(n * sizeof(double*));
    for (int row = 0; row < n; row++)
        result[row] = (double *) calloc(m , sizeof(double));
    return result;
 }
 
 
 double **createImplicitLaplacian(int dim, double r)
 {
    double **result = createZeroMatrix(dim, dim);
    for (int i = 0; i < dim; i++) {
        if (i == 0 || i == dim - 1) {
            result[i][i] = 1;
        } else {
            result[i][i] = 2.0*(1 + r);
            result[i][i-1] = -r;
            result[i][i+1] = -r;
        }
    }
    return result;
 }

 double *solveTriDiag(double **T, double *up, int n)
  {
    double a[n], b[n], c[n], c_prime[n], d[n], u_new[n];
    for (int i = 0; i < n; i++) {
        b[i] = T[i][i];
    }
    for (int i = 0; i < n-1; i++) {
        a[i+1] = T[i+1][i];
        c[i] = T[i][i+1];
    }
    c_prime[0] = c[0]/b[0];
    u_new[0] = (up[0])/b[0];
    for (int i = i; i <  n; i++) {
        c_prime[i] = c[i] / (b[i]  - c_prime[i-1] * a[i]);
        u_new[i] = ((up[i]) - u_new[i-1] * a[i])/(b[i] - c_prime[i-1] * a[i]);
    }

    for (int i = n-2; i >= 0; i--) {
        u_new[i] -= c_prime[i] * u_new[i+1];
    }

    for (int i = 0; i < n; i++) {
        up[i] = u_new[i];
    }    
    return up;   
}
 void freeMatrix(double **arr, int n, int m)
 {
    for (int row = 0; row < n; row++)
        free(arr[row]);
    free(arr);
 }
main(int argc, char * argv[]) {

  boil::timer.start();

  double xs[nx_s] = {}; 
  double xl[nx_l] = {};
  
  double T_l[nx_l][nt] = {{}};
  double T_s[nx_s][nt] = {{}};
  
  double t[nt] = {};

  double dx_l = (L_l - 0)/(nx_l - 1);
  double dx_s = (L_s - 0)/(nx_s - 1);
  double dt = (T - 0)/(nt - 1);
 
//////////////// ANALYTICAL SOLUTION //////////////////

  double Tf = 0;
  Tf = (q/lam_s + T_bc)/(1 - (lam_f/lam_s)) + T_bc;
  boil::oout << Tf << boil::endl;

  
  /*-----------------------+
  +	   DOMAIN 	   +
  +------------------------*/
  
  
  // (start) [:::::::::LIQUID:::::::::|/////SOLID/////]<<<Q (end) (x>+) 
  //         0                        x_in            L 
  
  for(int i = 0; i < nx_s; i++){
     xs[i] = i*dx_s  + x_in; 
  }
        
  for(int i = 0; i < nx_l; i++){
     xl[i] = i*dx_l;
  }
     
  for(int i = 0; i < nx_s; i++){
     T_s[i][0] = 273.0; 
  }
  for(int i = 0; i < nx_l; i++){
     T_l[i][0] = 273.0; 
  }
  
  double xs_min = xs[0];
  double xl_min = xl[0];
  double xs_max = xs[nx_s - 1];
  double xl_max = xl[nx_l - 1]; 

  double cp_s = rho_s*Cp_s;
  double cp_l = rho_l*Cp_f;
  
//  double del_li = (xl[nx_l - 2] - x_in);
//  double del_si = (xs[1] - x_in);
  
//  std::cout << del_si << boil::endl;
//  std::cout << del_li << boil::endl;
//  std::cout << dx_l << boil::endl;
//  std::cout << dx_s << boil::endl;
  
  /*-----------------------+
  +  BOUNDARY CONDITIONS   +			   
  +------------------------*/
  
  for(int i = 0; i < nt; i++){
  	T_l[0][i] = T_bc; //T = 100* C
  }
  
  //::::::::::: SOLVER ::::::::::://
 

  double Rl = dx_l/lam_f;
  double Rs = dx_s/lam_s;
  double kl = dt/(pow(dx_l,2));
  double ks = dt/(pow(dx_s,2));
    
  double r_s = (lam_s/cp_s)*ks;
  double r_l = (lam_f/cp_l)*kl;
    
  
  double b_s[nx_s] = {};
  double b_l[nx_l] = {};

  
//  boil::oout << T_l[0][0] << boil::endl;
  
  
  // initialisation of the new matrix
  
 /* 
  for(int i = 0; i < nx_s; i++){
      T_s_new[i][0] = T_s[i][0];
  }
  for(int i = 0; i < nx_l; i++){
      T_l_new[i][0] = T_l[i][0];
  }
  */
//  boil::oout << "new matrix initialised" << boil::endl;
  
  for(int j = 0; j < nt - 1; j++){
  boil::oout << "time loop started" << boil::endl;
  
  double **As = createImplicitLaplacian(nx_s, r_s);
  double **Al = createImplicitLaplacian(nx_l, r_l);
 
//  boil::oout << "Implict initialised" << boil::endl;

       
      for(int i = 1; i < nx_s - 1; i++){
       b_s[i] = 2.0*(1.0 - r_s) * T_s[i][j] + (r_s)*(T_s[i - 1][j]+ T_s[i + 1][j]);
      }
      b_s[0] = 1.0;
      b_s[nx_s -1] = 1.0;
      b_s[0] = r_s*(T_s[0][j]); //bc
      b_s[nx_s -1] = r_s*(T_s[nx_s - 1][j]); 
      
      double * bs = b_s;
//      double **bs1 = bs;
      
//      boil::oout << "d is initialised" << boil::endl;
  
      double * Ts = solveTriDiag(As, bs, nx_s);
//      boil::oout << "solid TDM solved, Crank Nicolson ran: " << j << " times" << boil::endl;

//    imposing boundary conditions
      T_s[nx_s-1][j+1] = (-q*dx_s)/lam_s + T_s[nx_s-2][j+1];
      
//    since one of our boundary condition is dependent on the temperature of the liquid side, we take it as

      T_s[0][j+1] = (Rs*T_l[nx_l-2][j] + Rl*T_s[1][j+1])/(Rl+Rs); // at the interface 
      T_l[nx_l - 1][j] = T_s[0][j+1];

      for(int i = 0; i < nx_s; i++){
//      boil::oout << Ts[i] << boil::endl;;
       }
       
      //updated values for the new time steps
//       boil::oout << "solid pointer to array done" << boil::endl;
      //fluid
      
         
 // at the interface;  
     

      for(int i = 1; i < nx_l - 1; i++){
        b_l[i] = 2.0*(1.0 - r_l) * T_l[i][j] + (r_l)*(T_l[i - 1][j]+ T_l[i + 1][j]);
      }
//      b_l[0] = 1.0;
//      b_l[nx_s -1] = 1.0;

       b_l[0] = r_l*(T_l[0][j]); //bc/
       b_l[nx_l -1] = r_l*(T_l[nx_l - 1][j]); //bc

;
      
//      boil::oout << "d is initialised for liquid" << boil::endl;
      double * bl = b_l;
//      double **bl1 = bl;
      
      double* Tl = solveTriDiag(Al, bl, nx_l);
      
      //updating the interaface temperature again
      T_s[0][j+1] = (Rs*T_l[nx_l-2][j+1] + Rl*T_s[1][j+1])/(Rl+Rs);
      T_l[nx_l - 1][j+1] = T_s[0][j+1];
      
      freeMatrix(As,nx_s,nx_s);
      freeMatrix(Al,nx_l,nx_l);

  }
  
  for(int i = 0; i < nt; i++){
 	std::cout << "liquid " << b_l[0]<< boil::endl;
 	std::cout << "solid " << b_s[0]<< boil::endl;
  }

  boil::oout << pi << boil::endl;
  boil::oout << Rl << boil::endl;
  boil::oout << Rs << boil::endl;
  boil::oout << r_l << boil::endl;
  boil::oout << r_s << boil::endl;
  
  boil::timer.stop();
  boil::timer.report();
}

