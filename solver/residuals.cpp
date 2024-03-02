#include "DG.h"

void u_quad(state_vec& uq, state_vec u[], int np, double phi[], int q, int k){
    uq.mass = 0, uq.mom_u=0, uq.mom_v=0, uq.energy=0;
    for(int j=0; j<np;j++){
        uq.mass = uq.mass + phi[q*np+j]*u[k*np+j].mass;
        uq.mom_u = uq.mom_u + phi[q*np+j]*u[k*np+j].mom_u;
        uq.mom_v = uq.mom_v + phi[q*np+j]*u[k*np+j].mom_v;
        uq.energy = uq.energy + phi[q*np+j]*u[k*np+j].energy;
    }
}

void element_residual(flux Residual[], double phi[], double dphi[], double J_inv[], state_vec u[],
                        double wq2[], double det[], int N_lin, int nq2, int np){
    //  N_elem is number of lineal elements
    state_vec u_q;
    double J1,J2,J3,J4;
    flux Fx, Fy;
    double dphi_x, dphi_y;
    //Loop over elements
    for(int k=0;k<N_lin;k++){
        // Loop to sum quad points
        J1 = J_inv[k*4]; J2 = J_inv[k*4+1]; J3 = J_inv[k*4+2]; J4 = J_inv[k*4+3];
        for(int q=0;q<nq2;q++){
            u_quad(u_q, u, np, phi, q, k); //compute u at quad point q for element k
            Fx = Flux_x(u_q);
            Fy = Flux_y(u_q);
            // Loop over unknowns of the element
            for(int i=0;i<np;i++){
                dphi_x = dphi[q*np*2+2*i], dphi_y = dphi[q*np*2+2*i+1];
                Residual[k*np+i].mass -=  (Fx.mass*(dphi_x*J1+dphi_y*J3) + Fy.mass*(dphi_x*J2+dphi_y*J4))*det[k]*wq2[q];
                Residual[k*np+i].mom_u -=  (Fx.mom_u*(dphi_x*J1+dphi_y*J3) + Fy.mom_u*(dphi_x*J2+dphi_y*J4))*det[k]*wq2[q];
                Residual[k*np+i].mom_v -=  (Fx.mom_v*(dphi_x*J1+dphi_y*J3) + Fy.mom_v*(dphi_x*J2+dphi_y*J4))*det[k]*wq2[q];
                Residual[k*np+i].energy -=  (Fx.energy*(dphi_x*J1+dphi_y*J3) + Fy.energy*(dphi_x*J2+dphi_y*J4))*det[k]*wq2[q];
            }
        }
    }
}

void edge_residual(flux Residual[], state_vec u[], int I2E[], double dl_int[], 
                    double In[], int N_int, int nq1, int np, double sq1[], double wq1[],
                     int B2E[], double dl_bound[], double Bn[], int N_bound, bool free_stream, 
                     double phi1[], double phi2[], double phi3[], int N_curved,
                     int N_lin){
    
    double phiR[nq1*np],phiL[nq1*np];
    flux flux_face;
    int index_L, index_R, local_indexL, local_indexR;
    int bgroup;
    double normal[2]; 
    double dl;
    state_vec uR, uL;
    double max_wsL[np], max_wsR[np];
    double J1,J2, deter;
    int curv_i;
    //Loop over the interior faces
    for(int f=0;f<N_int;f++){
        index_L = I2E[4*f];
        local_indexL = I2E[4*f+1];
        index_R = I2E[4*f+2];
        local_indexR = I2E[4*f+3];
        dl = dl_int[f];
        if(local_indexL==0){
            memcpy(phiL, phi1, sizeof(phiL));
        }
        else if(local_indexL==1){
            memcpy(phiL, phi2, sizeof(phiL));
        }
        else if(local_indexL==2){
            memcpy(phiL, phi3, sizeof(phiL));
        }
        if(local_indexR==0){
            memcpy(phiR, phi1, sizeof(phiR));
        }
        else if(local_indexR==1){
            memcpy(phiR, phi2, sizeof(phiR));
        }
        else if(local_indexR==2){
            memcpy(phiR, phi3, sizeof(phiR));
        }
        // print_matrix(phiL,0, np,np);
        for(int l=0;l<np;l++){
            max_wsL[l] = 0,max_wsR[l] = 0; 
        }
        normal[0] = In[2*f];
        normal[1] = In[2*f+1];
        //sum over the quad points
        for(int q=0;q<nq1;q++){
            u_quad(uL, u, np, phiL, q, index_L); //compute uL
            u_quad(uR, u, np, phiR, nq1-1-q, index_R); // compute uR
            flux_face = roe_flux(uR,uL, normal, true);
            for(int i=0;i<np;i++){
                //update Residual of the left element
                Residual[index_L*np+i].mass += phiL[q*np+i]*flux_face.mass*dl*wq1[q];
                Residual[index_L*np+i].mom_u += phiL[q*np+i]*flux_face.mom_u*dl*wq1[q];
                Residual[index_L*np+i].mom_v += phiL[q*np+i]*flux_face.mom_v*dl*wq1[q];
                Residual[index_L*np+i].energy += phiL[q*np+i]*flux_face.energy*dl*wq1[q];
                if(flux_face.wavespeed>max_wsL[i]){
                    max_wsL[i] = flux_face.wavespeed;
                }
                //update Residual of the right element
                Residual[index_R*np+i].mass -= phiR[(nq1-1-q)*np+i]*flux_face.mass*dl*wq1[q];
                Residual[index_R*np+i].mom_u -= phiR[(nq1-1-q)*np+i]*flux_face.mom_u*dl*wq1[q];
                Residual[index_R*np+i].mom_v -= phiR[(nq1-1-q)*np+i]*flux_face.mom_v*dl*wq1[q];
                Residual[index_R*np+i].energy -= phiR[(nq1-1-q)*np+i]*flux_face.energy*dl*wq1[q];
                if(flux_face.wavespeed>max_wsR[i]){
                    max_wsR[i] = flux_face.wavespeed;
                }
            }
        }
        for(int i=0;i<np;i++){
            Residual[index_L*np+i].wavespeed += max_wsL[i];
            Residual[index_R*np+i].wavespeed += max_wsR[i];
        }
    }
    //Loop over boundary faces
    for(int f=0;f<N_bound-N_curved;f++){
        index_L = B2E[3*f];
        local_indexL = B2E[3*f+1];
        bgroup = B2E[3*f+2];
        dl = dl_bound[f];
        normal[0] = Bn[2*f];
        normal[1] = Bn[2*f+1];
        if(local_indexL==0){
            memcpy(phiL, phi1, sizeof(phiL));
        }
        else if(local_indexL==1){
            memcpy(phiL, phi2, sizeof(phiL));
        }
        else if(local_indexL==2){
            memcpy(phiL, phi3, sizeof(phiL));
        }
        for(int l=0;l<np;l++){
            max_wsL[l] = 0; 
        }
        
        //sum over the quad points
        for(int q=0;q<nq1;q++){
            u_quad(uL, u, np, phiL, q, index_L); //compute uL
                if(free_stream){
                    full_state(flux_face, uL,normal, true);
                    
                }
                // Actually solve for the real boundary conditions
                else{
                    if(bgroup==0){
                        flux_face.mass = 0;
                        flux_face.mom_u = 0;
                        flux_face.mom_v = 0;
                        flux_face.energy = 0;
                        flux_face.wavespeed = 0;
                    }
                    else if(bgroup==1){
                        sub_outflow(flux_face,uL,normal, p_inf, true);
                    }
                    else if(bgroup==2){
                        inviscid(flux_face,uL,normal);
                    }
                    else if(bgroup==3){
                        inflow(flux_face,uL,normal,flow_angle,R,Tt,pt, true);
                    }
                }
            for(int i=0;i<np;i++){
                //update Residual of the left element
                Residual[index_L*np+i].mass += phiL[q*np+i]*flux_face.mass*dl*wq1[q];
                Residual[index_L*np+i].mom_u += phiL[q*np+i]*flux_face.mom_u*dl*wq1[q];
                Residual[index_L*np+i].mom_v += phiL[q*np+i]*flux_face.mom_v*dl*wq1[q];
                Residual[index_L*np+i].energy += phiL[q*np+i]*flux_face.energy*dl*wq1[q];
                if(flux_face.wavespeed>max_wsL[i]){
                    max_wsL[i] = flux_face.wavespeed;
                }
               }
        }
        for(int i=0;i<np;i++){
            Residual[index_L*np+i].wavespeed += max_wsL[i];
        }


    }
}

void element_residual_curved(flux Residual[], double phi[], double dphi[], double J_inv[], state_vec u[],
                        double wq2[], double det[],int N_lin, int N_curv, int nq2, int np){
    //  N_elem is number of curved elements
    state_vec u_q;
    double J1,J2,J3,J4;
    flux Fx, Fy;
    double dphi_x, dphi_y;
    //Loop over curved elements
    for(int k=0;k<N_curv;k++){
        // Loop to sum quad points
        for(int q=0;q<nq2;q++){
            u_quad(u_q, u, np, phi, q, k+N_lin); //compute u at quad point q for curved element k
            Fx = Flux_x(u_q);
            Fy = Flux_y(u_q);
            J1 = J_inv[k*nq2*4+4*q+0]; 
            J2 = J_inv[k*nq2*4+4*q+1]; 
            J3 = J_inv[k*nq2*4+4*q+2];
            J4 = J_inv[k*nq2*4+4*q+3];
            // Loop over unknowns of the element
            for(int i=0;i<np;i++){
                dphi_x = dphi[q*np*2+2*i], dphi_y = dphi[q*np*2+2*i+1];
                Residual[(N_lin+k)*np+i].mass -=  (Fx.mass*(dphi_x*J1+dphi_y*J3) + Fy.mass*(dphi_x*J2+dphi_y*J4))*det[k*nq2+q]*wq2[q];
                Residual[(N_lin+k)*np+i].mom_u -=  (Fx.mom_u*(dphi_x*J1+dphi_y*J3) + Fy.mom_u*(dphi_x*J2+dphi_y*J4))*det[k*nq2+q]*wq2[q];
                Residual[(N_lin+k)*np+i].mom_v -=  (Fx.mom_v*(dphi_x*J1+dphi_y*J3) + Fy.mom_v*(dphi_x*J2+dphi_y*J4))*det[k*nq2+q]*wq2[q];
                Residual[(N_lin+k)*np+i].energy -=  (Fx.energy*(dphi_x*J1+dphi_y*J3) + Fy.energy*(dphi_x*J2+dphi_y*J4))*det[k*nq2+q]*wq2[q];
            }
        }
    }
}


void edge_residual_curved(flux Residual[], state_vec u[], int N_elem,
                     int nq1, int np, double sq1[], double wq1[],
                     int B2E[], double dl_bound[],  int N_bound, bool free_stream, 
                     double phi_edge[], int N_curved, double Jc1[],double Jc2[],double Jc3[]){
    
    double phiL[np]; set_zero(phiL, np);
    flux flux_face;
    int index_L, local_indexL;
    int bgroup;
    int curv_i;
    double normal[2];
    double dl;
    double J1,J2;
    double deter;
    state_vec uL;
    double max_wsL[np];
    //Loop over boundary faces
    for(int f=N_bound-N_curved;f<N_bound;f++){
        index_L = B2E[3*f];
        local_indexL = B2E[3*f+1];
        bgroup = B2E[3*f+2];
        dl = dl_bound[f];
        curv_i = index_L-(N_elem-N_curved); // local index of the curved element
        for(int l=0;l<np;l++){
            max_wsL[l] = 0; 
        }
        //sum over the quad points
        for(int q=0;q<nq1;q++){
            if(local_indexL==0){
                J1 = Jc1[curv_i*nq1*2+2*q];
                J2 = Jc1[curv_i*nq1*2+2*q+1];
                if(np==1){
                    phiL[0] = phi_edge[q*np];
                }
                else if(np==3){
                    phiL[0]=0; phiL[1] = phi_edge[np*q]; phiL[2] = phi_edge[np*q+1];
                }
                else if(np==6){
                    phiL[0]=0; phiL[1] =0 ; phiL[2] = phi_edge[np*q];
                    phiL[3]=0; phiL[4] =phi_edge[np*q+1] ; phiL[5] = phi_edge[np*q+2];
                }
            }
            else if(local_indexL==1){
                J1 = Jc2[curv_i*nq1*2+2*q];
                J2 = Jc2[curv_i*nq1*2+2*q+1];
                if(np==1){
                    phiL[0] = phi_edge[q*np];
                }
                else if(np==3){
                    phiL[0]=phi_edge[np*q+1]; phiL[1] = 0; phiL[2] = phi_edge[np*q];
                }
                else if(np==6){
                    phiL[0]=phi_edge[np*q+2]; phiL[1] =0 ; phiL[2] = 0;
                    phiL[3]=phi_edge[np*q+1]; phiL[4] =0 ; phiL[5] = phi_edge[np*q];
                }
            }
            else if(local_indexL==2){
                J1 = Jc3[curv_i*nq1*2+2*q];
                J2 = Jc3[curv_i*nq1*2+2*q+1];
                if(np==1){
                    phiL[0] = phi_edge[q*np];
                }
                else if(np==3){
                    phiL[0]=phi_edge[np*q]; phiL[1] = phi_edge[np*q+1]; phiL[2] = 0;
                }
                else if(np==6){
                    phiL[0]=phi_edge[np*q]; phiL[1] = phi_edge[np*q+1] ; phiL[2] = phi_edge[np*q+2];
                    phiL[3]=0; phiL[4] =0 ; phiL[5] = 0;
                }
                
            }
            u_quad(uL, u, np, phiL, 0, index_L); //compute uL
            deter = sqrt(J1*J1+J2*J2);
            normal[0] = J2/deter;
            normal[1] = -J1/deter;
            if(free_stream){
                full_state(flux_face, uL,normal,true);
            }
            // Actually solve for the real boundary conditions
            else{
                    inviscid(flux_face,uL,normal);
            }
            for(int i=0;i<np;i++){
                //update Residual of the left element
                Residual[index_L*np+i].mass += phiL[i]*flux_face.mass*deter*wq1[q];
                Residual[index_L*np+i].mom_u += phiL[i]*flux_face.mom_u*deter*wq1[q];
                Residual[index_L*np+i].mom_v += phiL[i]*flux_face.mom_v*deter*wq1[q];
                Residual[index_L*np+i].energy += phiL[i]*flux_face.energy*deter*wq1[q];
                if(flux_face.wavespeed>max_wsL[i]){
                    max_wsL[i] = flux_face.wavespeed;
                }
            }
        }
        for(int i=0;i<np;i++){
            Residual[index_L*np+i].wavespeed += max_wsL[i];
        }
    }
}






