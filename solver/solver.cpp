
#include "DG.h"

void initialize(state_vec u[], int N_elem, int np){
    double rho,c,vx,vy,rhoE,T;
    T = 1.;
    rho = p_inf/(R*T);
    c = sqrt(gamm*R*T);
    vx =  M_inf*c;
    vy = 0.;
    rhoE = p_inf/(gamm-1) + 0.5*rho*vx*vx;
    for(int i=0;i<N_elem*np;i++){
        u[i].mass = rho;
        u[i].mom_u = rho*vx;
        u[i].mom_v = rho*vy;
        u[i].energy = rhoE;
    }
}

void reset_residual(flux Residual[],int N_elem, int np){
    for(int j=0; j<N_elem*np;j++){
        Residual[j].mass = 0;
        Residual[j].mom_u = 0;
        Residual[j].mom_v = 0;
        Residual[j].energy = 0;
        Residual[j].wavespeed = 0;
    }
}

void massinvmult(flux& F, double Mi[], flux Residual[], int l, int k, int init, int np, double det){
    F.mass = 0. ; F.mom_u = 0. ; F.mom_v = 0. ; F.energy = 0.;
    // print_matrix(Mi, init, np, np);
    for(int j=0;j<np;j++){
        F.mass = F.mass + Mi[init*np*np+np*l+j]*Residual[k*np+j].mass/det;
        F.mom_u = F.mom_u + Mi[init*np*np+np*l+j]*Residual[k*np+j].mom_u/det;
        F.mom_v = F.mom_v + Mi[init*np*np+np*l+j]*Residual[k*np+j].mom_v/det;
        F.energy = F.energy + Mi[init*np*np+np*l+j]*Residual[k*np+j].energy/det;
    }
}

void compute_Fn(flux Fn[], int N_elem, int  N_curv, int np, flux Residual[], double det[], double M_ref_inv[], double M_ki[]){
    flux F;
    int N_l = N_elem -N_curv;
    for(int j=0; j<N_elem*np;j++){
        Fn[j].mass = 0;
        Fn[j].mom_u = 0;
        Fn[j].mom_v = 0;
        Fn[j].energy = 0;
    }
    for(int k=0;k<N_l;k++){
        for(int l=0;l<np;l++){
            massinvmult(F, M_ref_inv, Residual, l, k ,0, np, det[k]);
            Fn[k*np+l].mass = -F.mass;
            Fn[k*np+l].mom_u = -F.mom_u;
            Fn[k*np+l].mom_v = -F.mom_v;
            Fn[k*np+l].energy = -F.energy;
        }
    }
    for(int k=0;k<N_curv;k++){
        for(int l=0;l<np;l++){
            massinvmult(F, M_ki, Residual, l, k+N_l, k, np, 1.);
            Fn[(k+N_l)*np+l].mass = -F.mass;
            Fn[(k+N_l)*np+l].mom_u = -F.mom_u;
            Fn[(k+N_l)*np+l].mom_v = -F.mom_v;
            Fn[(k+N_l)*np+l].energy = -F.energy;
        }
    }

}

void update_uF(state_vec uFE[], flux Fn[], double fac,double dt[], state_vec u[], int N_elem, int np){
    for(int j=0;j<N_elem;j++){
            for(int l=0;l<np;l++){
                uFE[j*np+l].mass = 0;
                uFE[j*np+l].mom_u = 0;
                uFE[j*np+l].mom_v = 0;
                uFE[j*np+l].energy = 0;
            }
    }
    for(int j=0;j<N_elem;j++){
            for(int l=0;l<np;l++){
                uFE[j*np+l].mass = u[j*np+l].mass + fac*dt[j*np+l]*Fn[j*np+l].mass;
                uFE[j*np+l].mom_u = u[j*np+l].mom_u + fac*dt[j*np+l]*Fn[j*np+l].mom_u;
                uFE[j*np+l].mom_v = u[j*np+l].mom_v + fac*dt[j*np+l]*Fn[j*np+l].mom_v;
                uFE[j*np+l].energy = u[j*np+l].energy + fac*dt[j*np+l]*Fn[j*np+l].energy;
            }
    }
}


double max_res(flux residual){
    double* max;
    double max_value;
    double values[4] = {abs(residual.mass),abs(residual.mom_u), abs(residual.mom_v),abs(residual.energy)};
    max = max_element(values,values+4);
    max_value = *max;
    return max_value;
}



// function to solve euler equation using DG
void solve_GD(int max_iter, bool free_stream, double tolerance,
                    int N_nodes, int N_elem, int N_curv, int N_int, int N_bound,
                    string files[6], int p, double CFL){
    string solution_file_name = "./Solutions/"+files[0]+"_sol_"+to_string(p)+".txt";
    string residual_file_name = "./Solutions/"+files[0]+"_res_"+to_string(p)+".txt";
    ofstream solution_file;
    solution_file.open(solution_file_name);
    ofstream residual_file;
    residual_file.open(residual_file_name);
    files[0] = "../curve_mesh/"+files[0]+".gri";
    int N_lin = N_elem - N_curv;
    double V[N_nodes*2]; //nodes coordinate
    int E[N_elem*3]; //element nodes
    int El[N_lin*3]; //element nodes
    int Ec[N_curv*10]; //element nodes
    int I2E[N_int*4]; // Interior faces to elements
    int B2E[N_bound*3]; // Boundary faces to elements
    double In[N_int*2]; // interior faces normal
    double Bn[N_bound*2]; // Boundary faces normal
    double Area[N_elem]; // Elements Area
    double P[N_elem]; set_zero(P,N_elem);
    double dl_int[N_int];
    double dl_bound[N_bound];
    double I_mid[2*N_int], B_mid[2*N_bound];
    // Assign values to the arrays
    read_file(V,El,Ec,I2E,B2E,In,Bn,P,Area,dl_int,dl_bound,files,N_int,N_bound,N_elem,N_curv,N_nodes,I_mid, B_mid);
    
    int np = (p+1)*(p+2)/2;   // # of coef per element
    // Quad points 1d and 2d
    int nq1; // # number of 1d points
    int nq1c; // # number of 1d points for curved elements
    if(p==0){
        nq1 = 1;
        nq1c = 3;
    }
    else if(p==1){
        nq1 = 2;
        nq1c = 4;
    }
    else if(p==2){
        nq1 = 3;
        nq1c = 5;
    }
    double sq1[nq1], wq1[nq1]; // nodes and weights 1d
    double sq1c[nq1c], wq1c[nq1c]; // nodes and weights 1d for curved elements
    quad1d(nq1,sq1,wq1);
    quad1d(nq1c,sq1c,wq1c);

    int nq2; // #number of 2d points
    int nq2c; // #number of 2d points for curved elements
    if(p==0){
        nq2 = 1;
        nq2c = 7;
    }
    else if(p==1){
        nq2 = 4;
        nq2c = 12;
    }
    else if(p==2){
        nq2 = 7;
        nq2c = 19;
    }
    double xyq[nq2*2], wq2[nq2]; // nodes and weights 2d
    double xyqc[nq2c*2], wq2c[nq2c]; // nodes and weights 2d for curved elements
    quad2d(nq2,xyq,wq2);
    quad2d(nq2c,xyqc,wq2c);

    double phi[nq2*np]; //phi at each quad point
    double dphi[nq2*np*2]; //dphi/dxi at each quad point
    for(int i=0; i<nq2;i++){
        basis_2d(np, phi, i, xyq[2*i], xyq[2*i+1]); // phi[0 to np] for the i quad point
        dphi_calc(np, dphi, i, xyq[2*i], xyq[2*i+1]); // same
    }
    double phic[nq2c*np]; //phi at each quad point for curved elements
    double dphic[nq2c*np*2]; //dphi/dxi at each quad point for curved elements
    
    for(int i=0; i<nq2c;i++){
        basis_2d(np, phic, i, xyqc[2*i], xyqc[2*i+1]); // phic[0 to np] for the i quad point
        dphi_calc(np, dphic, i, xyqc[2*i], xyqc[2*i+1]); // same
    }

    // phi at each edge
    double phi1[nq1*np], phi2[nq1*np], phi3[nq1*np]; //nq1 points per edge
    for(int q=0; q<nq1; q++){
        basis_2d(np, phi1, q, 1-sq1[q],sq1[q]);
        basis_2d(np, phi2, q, 0,1-sq1[q]);
        basis_2d(np, phi3, q, sq1[q],0);
    }

     // phi at each edge for curved elements
    double phi1c[nq1c*np], phi2c[nq1c*np], phi3c[nq1c*np]; //nq1 points per edge
    for(int q=0; q<nq1c; q++){
        basis_2d(np, phi1c, q, 1-sq1c[q],sq1c[q]);
        basis_2d(np, phi2c, q, 0,1-sq1c[q]);
        basis_2d(np, phi3c, q, sq1c[q],0);
    }

    // 1d phi at edge:
    double phic_1d[nq1c*np];
    for(int q=0; q<nq1c; q++){
        basis_1d(np, phic_1d, q, sq1c[q]);
    }
    

    // Jacobian and inverse lineal elements
    double J[N_lin*4]; //Jacobian
    double J_inv[N_lin*4]; // inverse of Jacobian
    double det[N_lin];
    Jacobian_linear(N_lin, El, V, J, J_inv, det); //compute the jacobian and its inverse per each element

    // Jacobian and inverse for curved elements at 2d quad points
    double Jc[N_curv*nq2c*4]; set_zero(Jc,N_curv*nq2c*4);//Jacobian
    double J_invc[N_curv*nq2c*4]; set_zero(J_invc,N_curv*nq2c*4);// inverse of Jacobian
    double detc[N_curv*nq2c];set_zero(detc,N_curv*nq2c);
    for(int i=0;i<N_curv;i++){
        for(int q=0;q<nq2c;q++){
            Jacobian_nonlinear(i, Ec, V, xyqc[2*q], xyqc[2*q+1], q, nq2c, Jc, J_invc, detc); //compute the jacobian and its inverse per each element
        }
    }
    

    // Jacobian at curved edges
    double Je1[N_curv*nq1c*2],dete1[N_curv*nq1c];
    set_zero(Je1,N_curv*nq1c*2); set_zero(dete1,N_curv*nq1c);
    double Je2[N_curv*nq1c*2],dete2[N_curv*nq1c];
    set_zero(Je2,N_curv*nq1c*2); set_zero(dete2,N_curv*nq1c);
    double Je3[N_curv*nq1c*2],dete3[N_curv*nq1c];
    set_zero(Je3,N_curv*nq1c*2); set_zero(dete3,N_curv*nq1c);
    for(int i=0;i<N_curv;i++){
        for(int q=0;q<nq1c;q++){
            Jacobian_edge(i, Ec, V, sq1c[q], q, nq1c, Je1, dete1, 0); //compute the jacobian and its inverse per each element
            Jacobian_edge(i, Ec, V, sq1c[q], q, nq1c, Je2, dete2, 1);
            Jacobian_edge(i, Ec, V, sq1c[q], q, nq1c, Je3, dete3, 2);
        }
    }
    // Mass matrix linear elements
    double M_ref[np*np]; 
    set_zero(M_ref,np*np);
    double M_ref_inv[np*np];
    set_zero(M_ref_inv,np*np);
    for(int i=0; i<np;i++){
        for(int j=0; j<np; j++){
            for(int q=0;q<nq2;q++){
                M_ref[i*np+j] +=  phi[q*np+i]*phi[q*np+j]*wq2[q];
            }
        }
    }
    // Compute the inverse of mass matrix and its inverse
    invert_matrix(M_ref, M_ref_inv, np,0);
    //Mass matrix curved elements 
    double M_curv[N_curv*np*np];set_zero(M_curv,N_curv*np*np);
    double M_curv_inv[N_curv*np*np];set_zero(M_curv_inv,N_curv*np*np);
    for(int k=0;k<N_curv;k++){
        M_k(np, nq2c, M_curv, k, phic, detc, wq2c);
        invert_matrix(M_curv, M_curv_inv, np, k);
    }
 
    
    // Define residual and some variables 
    flux Residual[N_elem*np];
    flux F0[N_elem*np];
    flux F1[N_elem*np];
    flux F2[N_elem*np];
    flux F3[N_elem*np];
    flux FF;
    double Mk[np*np];
    double dt[N_elem*np];
    double A_i,P_i,maxx,current;
    state_vec u[N_elem*np];
    state_vec uFE[N_elem*np];
    double local_max, global_max;
    /*
    //  ACTUALLY SOLVES FOR SPECIFIED MAX ITERATIONS
    */
   //  Initialize state vector u
    int contar=0;
    initialize(u, N_elem, np);
    cout << "++++++++++++START RUN+++++++++++" <<endl; 
    for(int t=0;t<max_iter;t++){

        // Compute F0 and time step
        reset_residual(Residual, N_elem, np);
        element_residual(Residual, phi, dphi, J_inv, u, wq2, det, N_lin, nq2,np);
        edge_residual(Residual, u, I2E, dl_int, In, N_int, nq1, np, sq1, wq1, B2E, dl_bound, Bn, N_bound, free_stream, phi1, phi2, phi3, N_curv, N_lin);
        // edge_residual_curved(Residual, u, N_elem, nq1c, np, sq1c, wq1c, B2E, dl_bound, N_bound, free_stream, phi1c, phi2c, phi3c, N_curv,Je1,Je2,Je3);
        edge_residual_curved(Residual, u, N_elem, nq1c, np, sq1c, wq1c, B2E, dl_bound, N_bound, free_stream, phic_1d, N_curv,Je1,Je2,Je3);
        // element_residual_curvedV2(Residual, phic, dphic,u, wq2c,N_lin, N_curv, nq2c, np, Ec, V, xyqc);
        element_residual_curved(Residual, phic, dphic,J_invc,u, wq2c, detc ,N_lin, N_curv, nq2c, np);
         // compute dt:
        for(int j=0;j<N_elem;j++){
            A_i = Area[j];
            P_i = P[j];
            maxx = 0;
            for(int q=0;q<np;q++){
                current = Residual[j*np+q].wavespeed;
                if(current>maxx){
                    maxx  = current;
                }
            }
            for(int l=0;l<np;l++){
                dt[j*np+l] = 2.*A_i*CFL/(maxx*P_i);
            }
        }

        compute_Fn(F0, N_elem, N_curv, np, Residual, det, M_ref_inv, M_curv_inv);
        update_uF(uFE, F0, 0.5,dt,u,N_elem,np);
        // Compute F1
        reset_residual(Residual, N_elem, np);
        element_residual(Residual, phi, dphi, J_inv, uFE, wq2, det, N_lin, nq2,np);
        edge_residual(Residual, uFE, I2E, dl_int, In, N_int, nq1, np, sq1, wq1, B2E, dl_bound, Bn, N_bound, free_stream, phi1, phi2, phi3, N_curv, N_lin);
        edge_residual_curved(Residual, uFE, N_elem, nq1c, np, sq1c, wq1c, B2E, dl_bound, N_bound, free_stream, phic_1d, N_curv,Je1,Je2,Je3);
        element_residual_curved(Residual, phic, dphic,J_invc,uFE, wq2c, detc ,N_lin, N_curv, nq2c, np);
        compute_Fn(F1, N_elem, N_curv, np, Residual, det, M_ref_inv, M_curv_inv);
        update_uF(uFE, F1, 0.5, dt,u,N_elem,np);

        // Compute F2
        reset_residual(Residual, N_elem, np);
        element_residual(Residual, phi, dphi, J_inv, uFE, wq2, det, N_lin, nq2,np);
        edge_residual(Residual, uFE, I2E, dl_int, In, N_int, nq1, np, sq1, wq1, B2E, dl_bound, Bn, N_bound, free_stream, phi1, phi2, phi3, N_curv, N_lin);
        edge_residual_curved(Residual, uFE, N_elem, nq1c, np, sq1c, wq1c, B2E, dl_bound, N_bound, free_stream, phic_1d, N_curv,Je1,Je2,Je3);
        element_residual_curved(Residual, phic, dphic,J_invc,uFE, wq2c, detc ,N_lin, N_curv, nq2c, np);
        compute_Fn(F2, N_elem, N_curv, np, Residual, det, M_ref_inv, M_curv_inv);
        update_uF(uFE, F2,1.,dt,u,N_elem,np);
        
        // Compute F3
        reset_residual(Residual, N_elem, np);
        element_residual(Residual, phi, dphi, J_inv, uFE, wq2, det, N_lin, nq2,np);
        edge_residual(Residual, uFE, I2E, dl_int, In, N_int, nq1, np, sq1, wq1, B2E, dl_bound, Bn, N_bound, free_stream, phi1, phi2, phi3, N_curv, N_lin);
        edge_residual_curved(Residual, uFE, N_elem, nq1c, np, sq1c, wq1c, B2E, dl_bound, N_bound, free_stream, phic_1d, N_curv,Je1,Je2,Je3);
        element_residual_curved(Residual, phic, dphic,J_invc,uFE, wq2c, detc ,N_lin, N_curv, nq2c, np);
        compute_Fn(F3, N_elem, N_curv, np, Residual, det, M_ref_inv, M_curv_inv);
        
        // compute u_n+1
        global_max = 0;
        for(int j=0;j<N_elem;j++){
            for(int l=0;l<np;l++){
                u[j*np+l].mass += (dt[j*np+l]/6.)*(F0[j*np+l].mass    +   2.*F1[j*np+l].mass     + 2.*F2[j*np+l].mass     +F3[j*np+l].mass);
                u[j*np+l].mom_u += (dt[j*np+l]/6.)*(F0[j*np+l].mom_u  +   2.*F1[j*np+l].mom_u    + 2.*F2[j*np+l].mom_u    +F3[j*np+l].mom_u);
                u[j*np+l].mom_v += (dt[j*np+l]/6.)*(F0[j*np+l].mom_v  +   2.*F1[j*np+l].mom_v    + 2.*F2[j*np+l].mom_v    +F3[j*np+l].mom_v);
                u[j*np+l].energy += (dt[j*np+l]/6.)*(F0[j*np+l].energy+   2.*F1[j*np+l].energy   + 2.*F2[j*np+l].energy   +F3[j*np+l].energy);
                local_max = max_res(Residual[j*np+l]);
                if(isnan(u[j].mass)){
                    cout << j << endl;
                    cout <<" SOME NAN" << endl;
                    exit(EXIT_FAILURE);
                }   
                if(local_max>global_max){
                    global_max = local_max;
                }
            }
        }
        
        cout << "Iteration " << t <<": Max Residual=" << global_max << endl;
        residual_file << global_max << '\n';
        if(global_max<tolerance){
            break;
        } 
    }

    // Write solution
    for(int i=0;i<N_elem*np;i++){
        solution_file << u[i].mass << '\t' << u[i].mom_u << '\t'
                        << u[i].mom_v << '\t' << u[i].energy << '\n';
    }
    residual_file.close();
    solution_file.close();
}

