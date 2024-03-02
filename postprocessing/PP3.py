import numpy as np
import matplotlib.pyplot as plt
import helper_fcns as hf


def plot_mesh(E,Ec,V,bot):
    x = np.zeros((4,2))
    plt.figure()
    for i in range(len(E)):
        x[0,0] = V[int(E[i,0]),0] ; x[0,1] = V[int(E[i,0]),1] 
        x[1,0] = V[int(E[i,1]),0] ; x[1,1] = V[int(E[i,1]),1] 
        x[2,0] = V[int(E[i,2]),0] ; x[2,1] = V[int(E[i,2]),1]
        x[3,0] = V[int(E[i,0]),0] ; x[3,1] = V[int(E[i,0]),1]
        plt.plot(x[:,0],x[:,1],'k-') 
    xm = np.zeros(2)    
    x = np.zeros((10,2))
    for i in range(len(Ec)):
        x[0,0] = V[int(Ec[i,0]),0] ; x[0,1] = V[int(Ec[i,0]),1] 
        x[1,0] = V[int(Ec[i,1]),0] ; x[1,1] = V[int(Ec[i,1]),1] 
        x[2,0] = V[int(Ec[i,2]),0] ; x[2,1] = V[int(Ec[i,2]),1]
        x[3,0] = V[int(Ec[i,3]),0] ; x[3,1] = V[int(Ec[i,3]),1]
        x[4,0] = V[int(Ec[i,6]),0] ; x[4,1] = V[int(Ec[i,6]),1]
        x[5,0] = V[int(Ec[i,8]),0] ; x[5,1] = V[int(Ec[i,8]),1]
        x[6,0] = V[int(Ec[i,9]),0] ; x[6,1] = V[int(Ec[i,9]),1]
        x[7,0] = V[int(Ec[i,7]),0] ; x[7,1] = V[int(Ec[i,7]),1]
        x[8,0] = V[int(Ec[i,4]),0] ; x[8,1] = V[int(Ec[i,4]),1]
        x[9,0] = V[int(Ec[i,0]),0] ; x[9,1] = V[int(Ec[i,0]),1]
        plt.plot(x[:,0],x[:,1],'k-')
        for j in  range(10):
            x[j,0] = V[int(Ec[i,j]),0]
            x[j,1] = V[int(Ec[i,j]),1]
        xc = (x[0,0]+x[3,0]+x[9,0])/3
        yc = (x[0,1]+x[3,1]+x[9,1])/3
        xm[0],xm[1] = hf.face_node(x,int(bot[i,1]),0.5)
        plt.plot(xm[0],xm[1],'ro')
        plt.text(xc,yc,str(i))
    plt.axis('equal')

def integrate_edge(V, local_face, sol,npp):
    nq = 5
    s, w = hf.quad1d(5)
    Cd = 0
    Cl = 0
    p_vec = np.zeros(5)
    x_vec = np.zeros(5)
    for q in range(5):
        J1,J2 = hf.Jacobian_edge(V, local_face, s[q])
        det = np.sqrt(J1*J1 + J2*J2)
        if(local_face==0):
            phi = hf.basis_2d(npp, 1-s[q],s[q])
        elif(local_face==1):
            phi = hf.basis_2d(npp, 0, 1-s[q])
        elif(local_face==2):
            phi = hf.basis_2d(npp, s[q], 0)
        nx = J2/det 
        ny = -J1/det
        rho = np.sum(phi*sol[:,0])
        mom_u = np.sum(phi*sol[:,1])
        mom_v = np.sum(phi*sol[:,2])
        E = np.sum(phi*sol[:,3])
        vx = mom_u/rho 
        vy = mom_v/rho
        vxb = vx - (vx*nx+vy*ny)*nx
        vyb = vy - (vx*nx+vy*ny)*ny
        p = (gamma-1)*(E-0.5*rho*(vx**2+vy**2))
        Cl += (p-p_inf)*ny*det*w[q]
        Cd += (p-p_inf)*nx*det*w[q]
        p_vec[q] = p
        x_vec[q],_ = hf.face_node(V,local_face, s[q]) 
    return Cl, Cd, p_vec, x_vec

def integrate_element(V, sol, npp, linear):
    rhoT = pt/(Tt*R)
    st = pt/(rhoT**gamma)
    xq,w = hf.quad2d()
    if(linear):
        J1 = V[1,0] - V[0,0]
        J2 = V[2,0] - V[0,0]
        J3 = V[1,1] - V[0,1]
        J4 = V[2,1] - V[0,1]
        det = J1*J4-J2*J3
        sq =0 ; integral=0
        for q in range(19):
            xi = xq[q,0] ; eta =xq[q,1]
            phi = hf.basis_2d(npp,xi,eta)
            rho = np.sum(phi*sol[:,0])
            mom_u = np.sum(phi*sol[:,1])
            mom_v = np.sum(phi*sol[:,2])
            E = np.sum(phi*sol[:,3])
            vx = mom_u/rho 
            vy = mom_v/rho 
            p = (gamma-1)*(E-0.5*rho*(vx**2+vy**2))
            sq = p/(rho**gamma)
            integral += (sq/st-1)**2*(det*w[q])
        A = 0
    else:
        sq =0 ; integral=0
        for q in range(19):
            xi = xq[q,0] ; eta =xq[q,1]
            det = hf.Jacobian_nonlinear(V,xi, eta)
            phi = hf.basis_2d(npp,xi,eta)
            rho = np.sum(phi*sol[:,0])
            mom_u = np.sum(phi*sol[:,1])
            mom_v = np.sum(phi*sol[:,2])
            E = np.sum(phi*sol[:,3])
            vx = mom_u/rho 
            vy = mom_v/rho 
            p = (gamma-1)*(E-0.5*rho*(vx**2+vy**2))
            sq = p/(rho**gamma)
            integral += (sq/st-1)**2*(det*w[q])
        A = hf.Area_curved(V,npp)
    return integral,A
        
def entropy(p, me):
    npp = int((p+1)*(p+2)/2)
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    sol = np.loadtxt('../solver/Solutions/bump'+str(me)+'_curv_sol_'+str(p)+'.txt')
    Area = np.loadtxt('../curve_mesh/Area'+str(me)+'.txt')
    Ec = mesh['Ec']
    El = mesh['El']
    V = mesh['V']
    init = len(El)
    A_i = 0 
    A_t = np.sum(Area[:init])
    E_t = 0
    E_i = 0
    for i in range(len(El)):
        x = np.zeros((3,2))
        for j in range(3):
            x[j,0] = V[El[i,j],0]; x[j,1] = V[El[i,j],1]
        sol_i = sol[i*npp:i*npp+npp,:]
        E_i,_ = integrate_element(x,sol_i,npp, True)
        E_t += E_i 
    for i in range(len(Ec)):
        x = np.zeros((10,2))
        for j in range(10):
            x[j,0] = V[Ec[i,j],0]; x[j,1] = V[Ec[i,j],1]
        sol_i = sol[(i+init)*npp:(i+init)*npp+npp,:]
        E_i, A_i = integrate_element(x,sol_i,npp, False)
        E_t += E_i
        A_t += A_i
    E_t = np.sqrt(E_t/A_t)
    return E_t

    
 
def coef(p, me):
    npp = int((p+1)*(p+2)/2)
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    B2E = np.loadtxt('../curve_mesh/B2E'+str(me)+'.txt')-1
    sol = np.loadtxt('../solver/Solutions/bump'+str(me)+'_curv_sol_'+str(p)+'.txt')
    Ec = mesh['Ec']
    El = mesh['El']
    dof = (len(Ec)+len(El))*npp
    V = mesh['V']
    bot = B2E[B2E[:,2]==0]
    init = len(El)
    bot[:,0] = bot[:,0] - init
    Cd =0
    Cl =0
    Cd_i =0
    Cl_i = 0
    p = np.zeros(5*len(Ec))
    x_bot = np.zeros(5*len(Ec))
    x = np.zeros((10,2))
    # Loop over edges
    for i in range(len(bot[:,0])):
        Elem_i = Ec[int(bot[i,0])]
        local_face = bot[i,1]
        for j in range(10):
            x[j,0] = V[Elem_i[j],0]; x[j,1] = V[Elem_i[j],1]
        sol_i = sol[(init+i)*npp:(init+i)*npp+npp,:]
        Cl_i, Cd_i,p[5*i:5*i+5],x_bot[5*i:5*i+5]  = integrate_edge(x, local_face, sol_i,npp)
        Cl += Cl_i 
        Cd += Cd_i
    den = gamma*p_inf*(M_inf**2)*h/2
    Cl = Cl/den
    Cd = Cd/den
    return dof,Cl, Cd,p, x_bot


def rate(x,y,name):
    slope, inter = np.polyfit(np.log(x),np.log(y),1)
    print("Convergence", name," :", slope)
    N = len(x)
    for i in range(N-1):
        pendiente = (np.log(y[i+1])-np.log(y[i]))/(np.log(x[i+1])-np.log(x[i]))
        print(pendiente)


def plot_stuff():
    mesh = [0]
    p = [0,1,2]
    Cd = np.zeros((3,3))
    Cl = np.zeros((3,3))
    E = np.zeros((3,3))
    Press = []
    x_bot = []
    dof = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            dof[i,j],Cl[i,j],Cd[i,j],p_l, x_i = coef(p[i], mesh[j])
            Press.append(p_l)
            x_bot.append(x_i)
            E[i,j] = entropy(p[i],mesh[j])
        rate(1/np.sqrt(dof[i,:]),np.abs(Cl[i,:]-Cl_teo),'Cl, p='+str(p[i]))
        rate(1/np.sqrt(dof[i,:]),np.abs(Cd[i,:]-Cd_teo),'Cd, p='+str(p[i]))
        rate(1/np.sqrt(dof[i,:]),np.abs(E[i,:]-E_teo),'E, p='+str(p[i]))

    col = ['orange', 'royalblue', 'g']
    plt.figure(figsize=(15,4))
    plt.subplot(1,3,1)
    for i in range(3):
        plt.semilogy(np.sqrt(dof[i,:]), np.abs(Cl[i,:]-Cl_teo),color=col[i], label="p="+str(p[i]))
    plt.legend()
    plt.grid(True)
    plt.xlabel(r'$\sqrt{dof}$')
    plt.ylabel(r'$C_l$  Error')
    plt.subplot(1,3,2)
    for i in range(3):
        plt.semilogy(np.sqrt(dof[i,:]), np.abs(Cd[i,:]-Cd_teo),color=col[i], label="p="+str(p[i]))
    plt.xlabel(r'$\sqrt{dof}$')
    plt.ylabel(r'$C_d$  Error')
    plt.grid(True)
    plt.legend()
    plt.subplot(1,3,3)
    for i in range(3):
        plt.semilogy(np.sqrt(dof[i,:]), np.abs(E[i,:]-E_teo),color=col[i], label="p="+str(p[i]))
    plt.xlabel(r'$\sqrt{dof}$')
    plt.ylabel(r'$E_s$  Error')
    plt.grid(True)
    plt.legend()
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('Error.png',bbox_inches='tight')

    plt.figure(figsize=(15,4))
    index_plot = [2,5,8]
    for i in range(3):
        plt.subplot(1,3,i+1)
        x = x_bot[index_plot[i]]
        P = Press[index_plot[i]]
        P = -(P-p_inf)/(gamma/2*p_inf*M_inf**2)
        index_sort = x.argsort()
        x = x[index_sort[::-1]]
        P = P[index_sort[::-1]]
        plt.plot(x,P,color = col[i],label='p='+str(p[i]))
        plt.xlabel(r'$x$')
        if(i==0):
            plt.ylabel(r'$-c_p$')
        plt.legend()
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('Pressure.png',bbox_inches='tight')

    






if __name__ == "__main__":
    gamma = 1.4
    R = 1
    p_inf = 1.
    M_inf = 0.5
    h = 0.0625
    Tt = 1+ (gamma-1)/2*M_inf**2 
    pt = Tt**(gamma/(gamma-1))
    Cl_teo = 1.537095
    Cd_teo = 2.94278e-6
    E_teo = 0
    
    #plot_stuff()
    dof0,Cl0,Cd0,p_l0, x_i0 = coef(0,0)
    dof1,Cl1,Cd1,p_l1, x_i1 = coef(1,0)
    dof2,Cl2,Cd2,p_l2, x_i2 = coef(2,0)

    plt.figure(figsize=(5,4))
    P = -(p_l0-p_inf)/(gamma/2*p_inf*M_inf**2)
    index_sort = x_i0.argsort()
    x = x_i0[index_sort[::-1]]
    P = P[index_sort[::-1]]
    plt.plot(x,P,label='p=0')

    P = -(p_l1-p_inf)/(gamma/2*p_inf*M_inf**2)
    index_sort = x_i1.argsort()
    x = x_i1[index_sort[::-1]]
    P = P[index_sort[::-1]]
    plt.plot(x,P,label='p=1')

    P = -(p_l2-p_inf)/(gamma/2*p_inf*M_inf**2)
    index_sort = x_i2.argsort()
    x = x_i2[index_sort[::-1]]
    P = P[index_sort[::-1]]
    plt.plot(x,P,label='p=2')

    plt.xlabel(r'$x$')
    plt.ylabel(r'$-c_p$')
    plt.legend()
    plt.subplots_adjust(hspace=0.25)
    plt.savefig('Pressure.png',bbox_inches='tight')
    #plt.show()
    