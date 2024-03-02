#include "DG.h"

//  q square
double q_2(state_vec u){
    double magnitude;
    magnitude = pow(u.mom_u/u.mass,2) + pow(u.mom_v/u.mass,2);
    return magnitude;
}

//  Total enthalpy per unit of mass
double H_calc(state_vec u){
    double Hval,vel_mag, p;
    vel_mag = q_2(u);
    p = (gamm-1)*(u.energy-0.5*u.mass*vel_mag);
    Hval = u.energy/u.mass + p/u.mass;
    return Hval;
}

// dot product between two vectors of length 2
double dot_prod(double v[2], double u[2]){
    double res;
    res = v[0]*u[0]+v[1]*u[1];
    return res;
}

// Flux x direction
flux Flux_x(state_vec u){
    flux F;
    double H,p;
    p = (gamm-1)*(u.energy-0.5*u.mass*q_2(u));
    H = H_calc(u);
    F.mass = u.mom_u;
    F.mom_u = pow(u.mom_u,2)/(u.mass)+p;
    F.mom_v = u.mom_u*u.mom_v/u.mass; 
    F.energy = u.mom_u*H;
    return F;
}

// Flux y direction
flux Flux_y(state_vec u){
    flux F;
    double H,p;
    p = (gamm-1)*(u.energy-0.5*u.mass*q_2(u));
    H = H_calc(u);
    F.mass = u.mom_v;
    F.mom_u = u.mom_u*u.mom_v/u.mass;
    F.mom_v = pow(u.mom_v,2)/(u.mass)+p;
    F.energy = u.mom_v*H;
    return F;
}

// compute the Roe flux given right and left vector states, and the normal.
flux roe_flux(state_vec uR, state_vec uL, double normal[2],bool flag_limiter=true){
    flux FL, FL_x, FL_y, FR, FR_x, FR_y, Roe;
    double v[2],H,HL,HR,q2,c,u;
    v[0] = (sqrt(uL.mass)*uL.mom_u/uL.mass+sqrt(uR.mass)*uR.mom_u/uR.mass)/(sqrt(uL.mass)+sqrt(uR.mass));
    v[1] = (sqrt(uL.mass)*uL.mom_v/uL.mass+sqrt(uR.mass)*uR.mom_v/uR.mass)/(sqrt(uL.mass)+sqrt(uR.mass));
    u = dot_prod(v,normal);
    HR = H_calc(uR);
    HL = H_calc(uL);
    H = (sqrt(uL.mass)*HL+sqrt(uR.mass)*HR)/(sqrt(uL.mass)+sqrt(uR.mass));
    q2 = pow(v[0],2)+pow(v[1],2);
    c = sqrt((gamm-1)*(H-0.5*q2));
    double lambda[4],s1,s2,G1,G2,C1,C2,Dmass, Dmom[2];
    Dmass = uR.mass-uL.mass ;
    lambda[0] = u+c;
    lambda[1] = u-c;
    lambda[2] = u;
    lambda[3] = u;
    if(flag_limiter){
        double eps=0.1*c;
        for(int i=0;i<4;i++){
            if(abs(lambda[i])<eps){
                lambda[i] = (pow(eps,2)+pow(lambda[i],2))/(2*eps);
            }
        }
    }

    s1 = 0.5*(abs(lambda[0])+abs(lambda[1]));
    s2 = 0.5*(abs(lambda[0])-abs(lambda[1]));
    Dmom[0] = uR.mom_u - uL.mom_u;
    Dmom[1] = uR.mom_v - uL.mom_v;
    G1 = (gamm-1)*(q2/2*Dmass-dot_prod(v,Dmom)+(uR.energy-uL.energy));
    G2 = -u*Dmass+dot_prod(Dmom, normal);
    C1 = G1/pow(c,2)*(s1-abs(lambda[2]))+G2/c*s2;
    C2 = G1/c*s2+(s1-abs(lambda[2]))*G2;
    FL_x = Flux_x(uL);
    FL_y = Flux_y(uL);
    FL.mass = FL_x.mass*normal[0] + FL_y.mass*normal[1];
    FL.mom_u = FL_x.mom_u*normal[0] + FL_y.mom_u*normal[1];
    FL.mom_v = FL_x.mom_v*normal[0] + FL_y.mom_v*normal[1];
    FL.energy = FL_x.energy*normal[0] + FL_y.energy*normal[1];     
    FR_x = Flux_x(uR);
    FR_y = Flux_y(uR);
    FR.mass = FR_x.mass*normal[0] + FR_y.mass*normal[1];
    FR.mom_u = FR_x.mom_u*normal[0] + FR_y.mom_u*normal[1];
    FR.mom_v = FR_x.mom_v*normal[0] + FR_y.mom_v*normal[1];
    FR.energy = FR_x.energy*normal[0] + FR_y.energy*normal[1];
    Roe.mass = 0.5*(FL.mass+FR.mass) - 0.5*(abs(lambda[2])*Dmass + C1);
    Roe.mom_u = 0.5*(FL.mom_u+FR.mom_u) - 0.5*(abs(lambda[2])*Dmom[0] + C1*v[0] + C2*normal[0]);
    Roe.mom_v = 0.5*(FL.mom_v+FR.mom_v) - 0.5*(abs(lambda[2])*Dmom[1] + C1*v[1] + C2*normal[1]);
    Roe.energy = 0.5*(FL.energy+FR.energy) - 0.5*(abs(lambda[2])*(uR.energy-uL.energy) + C1*H + C2*u);
    Roe.wavespeed = abs(u)+c;
    return Roe; 
}

// Inviscid Boundary condition
void inviscid(flux& flux_face,state_vec u, double normal[2]){
    double pb, vb[2],vx,vy,c;
    vx = u.mom_u/u.mass;
    vy = u.mom_v/u.mass;
    vb[0] = vx - (vx*normal[0]+vy*normal[1])*normal[0];
    vb[1] = vy - (vx*normal[0]+vy*normal[1])*normal[1];
    pb = (gamm-1)*(u.energy-0.5*u.mass*(vb[0]*vb[0]+vb[1]*vb[1]));
    c = sqrt(gamm*pb/u.mass);
    flux_face.mass = 0;
    flux_face.mom_u = pb*normal[0];
    flux_face.mom_v = pb*normal[1];
    flux_face.energy = 0;
    flux_face.wavespeed = abs(vb[0]*normal[0]+vb[1]*normal[1])+c;
}

// Inflow boundary condition
void inflow(flux& flux_face,state_vec u, double normal[2],
            double flow_angle, double R, double Tt, double pt, bool flag_limiter=true){
    double cp, up, Jp, dn, p, v_mag2;
    up = u.mom_u/u.mass*normal[0]+u.mom_v/u.mass*normal[1];
    v_mag2 = pow(u.mom_u/u.mass,2)+pow(u.mom_v/u.mass,2);
    p = (gamm-1)*(u.energy-0.5*u.mass*v_mag2);
    cp = sqrt(gamm*p/u.mass);
    Jp = up + 2*cp/(gamm-1);
    dn = normal[0]*cos(flow_angle) + normal[1]*sin(flow_angle);
    double a,b,c;
    a = pow(dn,2)*gamm*R*Tt-0.5*pow(Jp,2)*(gamm-1);
    b = 4.*dn*(gamm*R*Tt)/(gamm-1);
    c = 4.*gamm*R*Tt/(pow(gamm-1,2))-pow(Jp,2);
    double sol1,sol2,Mb;
    sol1 = (-b+sqrt(pow(b,2)-4*a*c))/(2*a);
    sol2 = (-b-sqrt(pow(b,2)-4*a*c))/(2*a);
    // assign value to Mb
    if(sol1>=0 && sol2>=0){
        if(sol1<sol2){
            Mb = sol1;
        }
        else{
            Mb = sol2;
        }
    }
    else if(sol1<0 && sol2>=0){
        Mb = sol2;
    }
    else if(sol1>=0 && sol2<0){
        Mb = sol1;
    }
    else{
        cout <<"Negative Mach number" << endl;
        cout << sol1<<endl;
        cout << sol2<<endl;
        cout <<pow(b,2)-4*a*c <<endl;
    }
    double Tb,pb, rhob, cb, vb[2], rhoEb;
    Tb = Tt/(1.+0.5*(gamm-1)*pow(Mb,2));
    pb = pt*pow(Tb/Tt,gamm/(gamm-1));
    rhob = pb/(R*Tb);
    cb = sqrt(gamm*pb/rhob);
    vb[0] = Mb*cb*cos(flow_angle);
    vb[1] = Mb*cb*sin(flow_angle);
    rhoEb = pb/(gamm-1)+0.5*rhob*(pow(vb[0],2)+pow(vb[1],2));
    //  Define state vector at boundary
    state_vec ub;
    ub.mass = rhob;
    ub.mom_u = rhob*vb[0];
    ub.mom_v = rhob*vb[1];
    ub.energy = rhoEb;
    flux_face = roe_flux(ub,ub,normal, flag_limiter);
}

// Outflow Boundary condition
void sub_outflow(flux& flux_face,state_vec u, double normal[2],
                double p_inf, bool flag_limiter=true){
    double cp, Jp, p, Sp,vp[2], v_mag2;
    vp[0] = u.mom_u/u.mass;
    vp[1] = u.mom_v/u.mass;
    v_mag2 = pow(vp[0],2)+pow(vp[1],2);
    p = (gamm-1)*(u.energy-0.5*u.mass*v_mag2);
    Sp = p/(pow(u.mass,gamm));
    cp = sqrt(gamm*p/u.mass);
    Jp = (vp[0]*normal[0]+vp[1]*normal[1])+2.*cp/(gamm-1);
    double rhob, pb, u_b, cb, vb[2], rhoEb;
    pb = p_inf;
    rhob = pow(pb/Sp,1./gamm);
    cb = sqrt(gamm*pb/rhob);
    u_b = Jp-2*cb/(gamm-1);
    vb[0] = vp[0]-(vp[0]*normal[0]+vp[1]*normal[1])*normal[0] + u_b*normal[0];
    vb[1] = vp[1]-(vp[0]*normal[0]+vp[1]*normal[1])*normal[1] + u_b*normal[1];
    rhoEb = pb/(gamm-1)+0.5*rhob*(pow(vb[0],2)+pow(vb[1],2));
    state_vec ub;
    ub.mass = rhob;
    ub.mom_u = rhob*vb[0];
    ub.mom_v = rhob*vb[1];
    ub.energy = rhoEb;
    flux_face = roe_flux(ub,ub,normal, flag_limiter);
}

void full_state(flux& flux_face,state_vec u, double normal[2],bool flag_limiter=true){
    double rho,c,vx,vy,rhoE,T;
    state_vec ub;
    T=1.;
    rho = p_inf/(R*T);
    c = sqrt(gamm*R*T);
    vx =  M_inf*c;
    vy = 0;
    rhoE = p_inf/(gamm-1) + 0.5*rho*vx*vx;
    ub.mass = rho;
    ub.mom_u = rho*vx;
    ub.mom_v = rho*vy;
    ub.energy = rhoE;
    flux_face = roe_flux(ub,u,normal, flag_limiter);
}


