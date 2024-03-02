import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
import scipy.interpolate
import helper_fcns as hf

def mach_0():
    me = 0
    p = 0
    npp =int((p+1)*(p+2)/2)
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    sol = np.loadtxt('../solver/Solutions/bump'+str(me)+'_curv_sol_'+str(p)+'.txt')
    Ec = mesh['Ec']
    El = mesh['El']
    V = mesh['V']
    init = len(El)
    nodes = []
    M_v = []
    triangles = []
    for i in range(len(El)):
        sol_i = sol[i*npp:i*npp+npp,:]
        rho = sol_i[:,0]
        mom_u = sol_i[:,1]
        mom_v = sol_i[:,2]
        rhoE = sol_i[:,3]
        vx = mom_u/rho
        vy = mom_v/rho
        p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
        c = np.sqrt(gamma*p/rho)
        mach_e = np.sqrt(vx**2+vy**2)/c
        n1 = mesh['El'][i,0]
        n2 = mesh['El'][i,1]
        n3 = mesh['El'][i,2]
        nT = [n1,n2,n3]
        triangles.append(nT)
        nodes.append(n1)
        nodes.append(n2)
        nodes.append(n3)
        M_v.append(mach_e)
        M_v.append(mach_e)
        M_v.append(mach_e)
    for i in range(len(Ec)):
        sol_i = sol[(i+init)*npp:(i+init)*npp+npp,:]
        rho = sol_i[:,0]
        mom_u = sol_i[:,1]
        mom_v = sol_i[:,2]
        rhoE = sol_i[:,3]
        vx = mom_u/rho
        vy = mom_v/rho
        p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
        c = np.sqrt(gamma*p/rho)
        mach_e = np.sqrt(vx**2+vy**2)/c
        n1 = mesh['Ec'][i,0]
        n2 = mesh['Ec'][i,3]
        n3 = mesh['Ec'][i,9]
        nT = [n1,n2,n3]
        triangles.append(nT)
        nodes.append(n1)
        nodes.append(n2)
        nodes.append(n3)
        M_v.append(mach_e)
        M_v.append(mach_e)
        M_v.append(mach_e)
    final_nodes = []
    Mach_ave = []
    nodes = np.asarray(nodes)
    M_v = np.asarray(M_v)
    while(True):
        if (len(nodes)==0):
            break
        n1 = nodes[0]
        index = np.where(nodes==n1)
        final_nodes.append(n1)
        Mach_ave.append(np.average(M_v[index]))
        nodes = np.delete(nodes,index)
        M_v = np.delete(M_v,index)
    x_v = np.zeros(len(final_nodes))
    y_v = np.zeros(len(final_nodes))
    Mach_ave = np.asarray(Mach_ave)
    final_nodes = np.asarray(final_nodes)
    Mach = np.zeros(len(final_nodes))
    for i in range(len(final_nodes)):
        x_v[int(final_nodes[i])] = V[int(final_nodes[i]),0];y_v[int(final_nodes[i])] = V[int(final_nodes[i]),1]
        Mach[int(final_nodes[i])] = Mach_ave[i]

    
    triang = mtri.Triangulation(x_v,y_v, triangles)
    fig, axes = plt.subplots(figsize=(10,4))
    axes.set_ylabel('y')
    axes.set_xticklabels([])
    axes.axis('equal')
    im =axes.tricontourf(triang, Mach, 15, cmap=cm.coolwarm,vmin=0.425,vmax=0.8)
    plt.colorbar(im,ax=axes, orientation='horizontal')
    axes.triplot(triang,'k-',lw=0.1)
    plt.ylim([0,0.8])
    plt.axis('off')
    plt.savefig('Mach0.png',bbox_inches='tight')


def mach_1():
    me = 0
    p = 1
    npp =int((p+1)*(p+2)/2)
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    sol = np.loadtxt('../solver/Solutions/bump'+str(me)+'_curv_sol_'+str(p)+'.txt')
    Ec = mesh['Ec']
    El = mesh['El']
    V = mesh['V']
    init = len(El)
    nodes = []
    M_v = []
    triangles = []
    for i in range(len(El)):
        sol_i = sol[i*npp:i*npp+npp,:]
        rho = sol_i[:,0]
        mom_u = sol_i[:,1]
        mom_v = sol_i[:,2]
        rhoE = sol_i[:,3]
        vx = mom_u/rho
        vy = mom_v/rho
        p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
        c = np.sqrt(gamma*p/rho)
        mach_e = np.sqrt(vx**2+vy**2)/c
        n1 = mesh['El'][i,0]
        n2 = mesh['El'][i,1]
        n3 = mesh['El'][i,2]
        nT = [n1,n2,n3]
        triangles.append(nT)
        nodes.append(n1)
        nodes.append(n2)
        nodes.append(n3)
        M_v.append(mach_e[0])
        M_v.append(mach_e[1])
        M_v.append(mach_e[2])
    for i in range(len(Ec)):
        sol_i = sol[(i+init)*npp:(i+init)*npp+npp,:]
        rho = sol_i[:,0]
        mom_u = sol_i[:,1]
        mom_v = sol_i[:,2]
        rhoE = sol_i[:,3]
        vx = mom_u/rho
        vy = mom_v/rho
        p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
        c = np.sqrt(gamma*p/rho)
        mach_e = np.sqrt(vx**2+vy**2)/c
        n1 = mesh['Ec'][i,0]
        n2 = mesh['Ec'][i,3]
        n3 = mesh['Ec'][i,9]
        nT = [n1,n2,n3]
        triangles.append(nT)
        nodes.append(n1)
        nodes.append(n2)
        nodes.append(n3)
        M_v.append(mach_e[0])
        M_v.append(mach_e[1])
        M_v.append(mach_e[2])
    final_nodes = []
    Mach_ave = []
    nodes = np.asarray(nodes)
    M_v = np.asarray(M_v)
    while(True):
        if (len(nodes)==0):
            break
        n1 = nodes[0]
        index = np.where(nodes==n1)
        final_nodes.append(n1)
        Mach_ave.append(np.average(M_v[index]))
        nodes = np.delete(nodes,index)
        M_v = np.delete(M_v,index)
    x_v = np.zeros(len(final_nodes))
    y_v = np.zeros(len(final_nodes))
    Mach_ave = np.asarray(Mach_ave)
    final_nodes = np.asarray(final_nodes)
    Mach = np.zeros(len(final_nodes))
    for i in range(len(final_nodes)):
        x_v[int(final_nodes[i])] = V[int(final_nodes[i]),0];y_v[int(final_nodes[i])] = V[int(final_nodes[i]),1]
        Mach[int(final_nodes[i])] = Mach_ave[i]
    
    triang = mtri.Triangulation(x_v,y_v, triangles)
    fig, axes = plt.subplots(figsize=(10,4))
    axes.set_ylabel('y')
    axes.set_xticklabels([])
    axes.axis('equal')
    im =axes.tricontourf(triang, Mach, 15, cmap=cm.coolwarm)
    # axes.triplot(triang,'k-',lw=0.1)
    plt.ylim([0,0.8])
    plt.axis('off')
    plt.colorbar(im, orientation='horizontal')
    plt.savefig('Mach1.png',bbox_inches='tight')
    return triang

def mach_2():
    me = 0
    p = 2
    npp =int((p+1)*(p+2)/2)
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    sol = np.loadtxt('../solver/Solutions/bump'+str(me)+'_curv_sol_'+str(p)+'.txt')
    B2E = np.loadtxt('../curve_mesh/B2E'+str(me)+'.txt')-1
    I2E = np.loadtxt('../curve_mesh/I2E'+str(me)+'.txt')-1
    final_node = get_linearnodes()
    Ec = mesh['Ec']
    El = mesh['El']
    V = mesh['V']
    linear_n = get_linearnodes()
    c_n = linear_n
    init = len(El)
    clock = [2,0,1]
    anti = [1,2,0]
    x =[]
    y = []
    nodes =[]
    Mach = []
    for i in range(len(I2E)):
        E_left = I2E[i,0]
        l_f = I2E[i,1]
        if(E_left<init):
            E_i = El[int(E_left)]
            n1 = E_i[anti[int(l_f)]]
            n2 = E_i[clock[int(l_f)]]
            n3 = c_n ;c_n =+1
            x1 = V[n1,0]; y1 = V[n1,1]
            x2 = V[n2,0]; y2 = V[n2,1]
            xm = (x1+x2)/2 ; ym= (y1+y2)/2
            nodes.append(n1);nodes.append(n2); nodes.append(n3)
            x.append(x1),x.append(x2);x.append(xm)
            y.append(y1);y.append(y2);y.append(ym)
            sol_i = sol[int(E_left)*npp:int(E_left)*npp+npp,:]
            rho = sol_i[:,0]
            mom_u = sol_i[:,1]
            mom_v = sol_i[:,2]
            rhoE = sol_i[:,3]
            vx = mom_u/rho
            vy = mom_v/rho
            p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
            c = np.sqrt(gamma*p/rho)
            mach_e = np.sqrt(vx**2+vy**2)/c
            if(l_f==0):
                Mach.append(mach_e[2])
                Mach.append(mach_e[5])
                Mach.append(mach_e[4])
            elif(l_f==1):
                Mach.append(mach_e[5])
                Mach.append(mach_e[0])
                Mach.append(mach_e[3])
            else:
                Mach.append(mach_e[0])
                Mach.append(mach_e[2])
                Mach.append(mach_e[1])

        else:
            E_i = Ec[int(E_left-init)]
            columns = [1,2,4,5,6,7,8]
            E_i = np.delete(E_i,columns)
            n1 = E_i[anti[int(l_f)]]
            n2 = E_i[clock[int(l_f)]]
            n3 = c_n ;c_n =+1
            x1 = V[n1,0]; y1 = V[n1,1]
            x2 = V[n2,0]; y2 = V[n2,1]
            xm = (x1+x2)/2 ; ym= (y1+y2)/2
            nodes.append(n1);nodes.append(n2); nodes.append(n3)
            x.append(x1),x.append(x2);x.append(xm)
            y.append(y1);y.append(y2);y.append(ym)
            sol_i = sol[int(E_left)*npp:int(E_left)*npp+npp,:]
            rho = sol_i[:,0]
            mom_u = sol_i[:,1]
            mom_v = sol_i[:,2]
            rhoE = sol_i[:,3]
            vx = mom_u/rho
            vy = mom_v/rho
            p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
            c = np.sqrt(gamma*p/rho)
            mach_e = np.sqrt(vx**2+vy**2)/c
            if(l_f==0):
                Mach.append(mach_e[2])
                Mach.append(mach_e[5])
                Mach.append(mach_e[4])
            elif(l_f==1):
                Mach.append(mach_e[5])
                Mach.append(mach_e[0])
                Mach.append(mach_e[3])
            else:
                Mach.append(mach_e[0])
                Mach.append(mach_e[2])
                Mach.append(mach_e[1])
    for i in range(len(B2E)):
        E_left = B2E[i,0]
        l_f = B2E[i,1]
        bound = B2E[i,2]
        if(bound==0):
            x_c = np.zeros((10,2))
            E_i = Ec[int(E_left-init)]
            for i in range(10):
                x_c[i,0] = V[int(E_i[i]),0];x_c[i,1] = V[int(E_i[i]),1]
            xm, ym = hf.face_node(x_c,l_f,0.5)
            columns = [1,2,4,5,6,7,8]
            E_i = np.delete(E_i,columns)
            n1 = E_i[anti[int(l_f)]]
            n2 = E_i[clock[int(l_f)]]
            n3 = c_n ;c_n =+1
            x1 = V[n1,0]; y1 = V[n1,1]
            x2 = V[n2,0]; y2 = V[n2,1]
            nodes.append(n1);nodes.append(n2); nodes.append(n3)
            x.append(x1),x.append(x2);x.append(xm)
            y.append(y1);y.append(y2);y.append(ym)
            sol_i = sol[int(E_left)*npp:int(E_left)*npp+npp,:]
            rho = sol_i[:,0]
            mom_u = sol_i[:,1]
            mom_v = sol_i[:,2]
            rhoE = sol_i[:,3]
            vx = mom_u/rho
            vy = mom_v/rho
            p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
            c = np.sqrt(gamma*p/rho)
            mach_e = np.sqrt(vx**2+vy**2)/c
            if(l_f==0):
                Mach.append(mach_e[2])
                Mach.append(mach_e[5])
                Mach.append(mach_e[4])
            elif(l_f==1):
                Mach.append(mach_e[5])
                Mach.append(mach_e[0])
                Mach.append(mach_e[3])
            else:
                Mach.append(mach_e[0])
                Mach.append(mach_e[2])
                Mach.append(mach_e[1])
        else:
            E_i = El[int(E_left)]
            n1 = E_i[anti[int(l_f)]]
            n2 = E_i[clock[int(l_f)]]
            n3 = c_n ;c_n =+1
            x1 = V[n1,0]; y1 = V[n1,1]
            x2 = V[n2,0]; y2 = V[n2,1]
            xm = (x1+x2)/2 ; ym= (y1+y2)/2
            nodes.append(n1);nodes.append(n2); nodes.append(n3)
            x.append(x1),x.append(x2);x.append(xm)
            y.append(y1);y.append(y2);y.append(ym)
            sol_i = sol[int(E_left)*npp:int(E_left)*npp+npp,:]
            rho = sol_i[:,0]
            mom_u = sol_i[:,1]
            mom_v = sol_i[:,2]
            rhoE = sol_i[:,3]
            vx = mom_u/rho
            vy = mom_v/rho
            p = (gamma-1)*(rhoE-0.5*rho*(vx**2+vy**2))
            c = np.sqrt(gamma*p/rho)
            mach_e = np.sqrt(vx**2+vy**2)/c
            if(l_f==0):
                Mach.append(mach_e[2])
                Mach.append(mach_e[5])
                Mach.append(mach_e[4])
            elif(l_f==1):
                Mach.append(mach_e[5])
                Mach.append(mach_e[0])
                Mach.append(mach_e[3])
            else:
                Mach.append(mach_e[0])
                Mach.append(mach_e[2])
                Mach.append(mach_e[1])

    x_b = []
    y_b = []
    tol =0.002
    for i in range(len(x)):
        if(np.abs(y[i]- 0.0625*np.exp(-25.*x[i]**2))<tol):
            x_b.append(x[i])
            y_b.append(y[i])
    x_b = np.asarray(x_b)
    y_b = np.asarray(y_b)
    index = x_b.argsort()
    x_b = x_b[index]
    y_b = y_b[index]
    #triang2 = mach_1()
    nodes = np.asarray(nodes)
    Mach = np.asarray(Mach)
    x = np.asarray(x)
    y = np.asarray(y)
    triang = mtri.Triangulation(x,y)
    fig, axes = plt.subplots(figsize=(10,4))
    axes.set_ylabel('y')
    axes.set_xticklabels([])
    axes.axis('equal')
    im =axes.tricontourf(triang, Mach, 15, cmap=cm.coolwarm,vmin=0.425,vmax=0.8)
    axes.triplot(triang,'k-',lw=0.1)
    for i in range(len(x_b)-1):
        x_f = np.linspace(x_b[i],x_b[i+1],10)
        y_f = np.linspace(y_b[i],y_b[i+1],10)
        plt.fill_between(x_f,y_f,color ='white')
    plt.ylim([0,0.8])
    #plt.clim(0.425,0.775)
    plt.axis('off')
    plt.colorbar(im, ax=axes,orientation='horizontal')
    plt.savefig('Mach2.png',bbox_inches='tight')
    

   

def get_linearnodes():
    me = 0
    mesh = hf.readgri('../curve_mesh/bump'+str(me)+'_curv.gri')
    Ec = mesh['Ec']
    El = mesh['El']
    nodes = [] 
    for i in range(len(El)):
        for j in range(3):
            nodes.append(El[i,j])
    for i in range(len(Ec)):
        nodes.append(Ec[i,0])
        nodes.append(Ec[i,3])
        nodes.append(Ec[i,9])
    nodes = np.asarray(nodes)
    return np.amax(nodes)

if __name__ == "__main__":
    gamma = 1.4
    R = 1
    p_inf = 1.
    M_inf = 0.5
    h = 0.0625
    Tt = 1+ (gamma-1)/2*M_inf**2 
    pt = Tt**(gamma/(gamma-1))
    
    mach_0()
    mach_2()
    #plt.show()