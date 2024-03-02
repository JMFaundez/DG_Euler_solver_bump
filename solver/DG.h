#ifndef MYHEADER_H
#define MYHEADER_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <armadillo> // g++ example.cpp -o example -O2 -larmadillo

using namespace std;
using namespace arma;

// Flow constants
const double gamm = 1.4;
const double M_inf = 0.5;
const double R = 1.0;
const double p_inf = 1.0;
const double flow_angle = 0;
const double Tt = 1.0+(gamm-1)/2.*pow(M_inf,2);
const double pt = pow(Tt,gamm/(gamm-1));


// state vector u [ rho, rho*u, rho*v, rho*E]
struct state_vec{
    double mass;
    double mom_u;
    double mom_v;
    double energy;
};

// Flux of the conserved quantities  
struct flux{
    double mass;
    double mom_u;
    double mom_v;
    double energy;
    double wavespeed;
};

// Read grid
void read_file(double V[], int El[], int Ec[],int I2E[], int B2E[],
                double In[],double Bn[],double P[],double Area[],
                double dl_int[], double dl_bound[],
                string namefile[6],int N_int, int N_bound,
                int N_elem, int N_curv, int N_nodes, double I_mid[], double B_mid[]);

//Matrix Operations
void print_matrix(double Ma[], int start, int n_row, int n_col);
void print_matrix_int(int Ma[], int start, int n_row, int n_col);
void invert_matrix(double M[], double M_inv[], int n_row, int init);
void set_zero(double M[], int elem);



//function to compute fluxes
flux Flux_x(state_vec u);
flux Flux_y(state_vec u);
flux roe_flux(state_vec uR, state_vec uL, double normal[2],bool flag_limiter);

// Boundary conditions
void inviscid(flux& flux_face,state_vec u, double normal[2]);
void inflow(flux& flux_face,state_vec u, double normal[2],double flow_angle, double R, double Tt, double pt, bool flag_limiter);
void sub_outflow(flux& flux_face,state_vec u, double normal[2],double p_inf, bool flag_limiter);
void full_state(flux& flux_face,state_vec u, double normal[2], bool flag_limiter);


//Residual Functions
void element_residual(flux Residual[], double phi[], double dphi[], double J_inv[], state_vec u[],
                        double wq2[], double det[], int N_lin, int nq2, int np);
void element_residual_curved(flux Residual[], double phi[], double dphi[], double J_inv[], state_vec u[],
                        double wq2[], double det[],int N_lin, int N_elem, int nq2, int np);
void edge_residual(flux Residual[], state_vec u[], int I2E[], double dl_int[], 
                    double In[], int N_int, int nq1, int np, double sq1[], double wq1[],
                     int B2E[], double dl_bound[], double Bn[], int N_bound, bool free_stream, 
                     double phi1[], double phi2[], double phi3[], int N_curved, int N_lin);
void edge_residual_curved(flux Residual[], state_vec u[], int N_elem,
                     int nq1, int np, double sq1[], double wq1[],
                     int B2E[], double dl_bound[],  int N_bound, bool free_stream, 
                     double phi_edge[], int N_curved, double Jc1[],double Jc2[],double Jc3[]);

// Solver functions 
void solve_GD(int max_iter, bool free_stream, double tolerance,int N_nodes, int N_elem, int N_curv, int N_int, int N_bound, string files[6], int p, double CFL);
double max_res(flux residual);
void update_uF(state_vec uFE[], flux Fn[], double fac, double dt[], state_vec u[], int N_elem, int np);
void compute_Fn(flux Fn[], int N_elem, int N_curv, int np, flux Residual[], double det[], double M_ref_inv[], double M_k);
void massinvmult(flux& F, double Mi[], flux Residual[], int i, int k,int init, int np, double det);
void reset_residual(flux Residual[],int N_elem, int np);
void initialize(state_vec u[], int N_elem, int np);


// Quad points, Jacobian, basis functions and Mass matrix
void quad2d(int nq, double x[], double w[]);
void quad1d(int q, double x[], double w[]);
void basis_2d(int np, double phi[],int quad_i, double xi, double eta);
void dphi_calc(int np, double dphi[],int quad_i, double xi, double eta);
void basis_1d(int np, double phi[],int quad_i, double s);
void Jacobian_linear(int N_elem, int E[], double V[], double J[], double J_inv[], double det[]);
void basis_geom_2d(double phi[], double dphi[], double xi, double eta);
void basis_geom_1d(double phi[], double dphi[], double s);
void Jacobian_nonlinear(int elem_i, int E[], double V[], double xi, double eta, int q,  int nq, double J[], double J_inv[], double det[]);
void Jacobian_edge(int elem_i, int E[], double V[], double s, int q,  int nq, double J[], double det[], int local);
void M_k(int np, int nq, double M[], int elem_i, double phi[],double det[], double w[]);
#endif /* MYHEADER_H */
