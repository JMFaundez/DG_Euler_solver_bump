import numpy as np

def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    V = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    B = []; Bname = []
    for i in range(NB):
        s = f.readline().split(); Nb = int(s[0]); Bname.append(s[2])
        Bi = np.array([[int(s)-1 for s in f.readline().split()] for n in range(Nb)])
        B.append(Bi)
    # read linear elements
    Ne0 = 0; El = []
    s = f.readline().split(); ne = int(s[0])
    Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
    El = Ei if (Ne0==0) else np.concatenate((El,Ei), axis=0)
    # read curved elements
    Ne0 = 0; Ec = []
    s = f.readline().split(); ne = int(s[0])
    Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
    Ec = Ei if (Ne0==0) else np.concatenate((Ec,Ei), axis=0)
    f.close()
    Mesh = {'V':V, 'El':El, 'Ec': Ec, 'B':B, 'Bname':Bname }
    return Mesh

def basis_2d(npp,xi, eta):
    phi = np.zeros(npp)
    if(npp==1):
        phi[0] = 1.
    elif(npp == 3):
        phi[0] = 1. - xi - eta
        phi[1] = xi
        phi[2] = eta
    elif(npp == 6):
        phi[0] = 1 - 3*xi + 2*(xi**2) - 3*eta + 4*xi*eta + 2*(eta**2)
        phi[1] = 4*xi - 4*(xi**2)-4*xi*eta
        phi[2] = -xi+2*(xi**2)
        phi[3] = 4*eta - 4*xi*eta-4*(eta**2)
        phi[4] = 4*xi*eta
        phi[5] = -eta + 2*(eta**2)
    return phi    

def dbasis_2d(npp, xi, eta):
    dphi = np.zeros(npp*2)
    if(npp==1):
        dphi[0] = 0; dphi[1] = 0
    elif(npp==3):
        dphi[0] = -1. ; dphi[1] = -1.
        dphi[2] = 1.; dphi[3] = 0.
        dphi[4] = 0.; dphi[5] = 1.
    elif(npp==6):
        dphi[0] = -3.+4.*xi+4.*eta ; dphi[1] = -3.+4.*xi+4.*eta
        dphi[2] = 4. - 8.*xi - 4*eta; dphi[3] = -4.*xi
        dphi[4] = -1.+4*xi; dphi[5] = 0.
        dphi[6] = -4.*eta ; dphi[7] = 4.-4.*xi-8.*eta
        dphi[8] = 4.*eta; dphi[9] = 4.*xi
        dphi[10] = 0; dphi[11] = -1.+4*eta
    return dphi

def basis_geom_1d(s):
    phi = np.zeros(4)
    dphi = np.zeros(4)
    phi[0] = 1.5*(-3.0*s + 1.0)*(s - 1.)*(s - 0.666666666666667)
    phi[1] = 13.5*s*(s - 1.)*(s - 0.666666666666667)
    phi[2] = -13.5*s*(s - 1.)*(s - 0.333333333333333)
    phi[3] = 4.5*s*(s - 0.666666666666667)*(s - 0.333333333333333)
    dphi[0] = (-4.5*s + 1.5)*(s - 0.666666666666667) + (-4.5*s + 4.5)*(s - 0.666666666666667) + (-3.0*s + 1.0)*(1.5*s - 1.5)
    dphi[1] = 13.5*s*(s - 1.) + 13.5*s*(s - 0.666666666666667) + (s - 0.666666666666667)*(13.5*s - 13.5)
    dphi[2] = -13.5*s*(s - 1.) - 13.5*s*(s - 0.333333333333333) + (-13.5*s + 13.5)*(s - 0.333333333333333)
    dphi[3] = 4.5*s*(s - 0.666666666666667) + 4.5*s*(s - 0.333333333333333) + (s - 0.666666666666667)*(4.5*s - 1.5)
    return  phi, dphi


def basis_geom_2d(xi, eta):
    phi = np.zeros(10)
    dphi = np.zeros((10,2))
    coef =  np.array([
                            [1.,    -5.5,   9.,     -4.5,   -5.5,   18.,    -13.5,  9.,     -13.5,  -4.5],
                            [0,     9.,     -22.5,  13.5,   0,      -22.5,  27,     0,      13.5,   0 ],
                            [0,     -4.5,   18,     -13.5,  0,      4.5,    -13.5,  0,      0,      0],
                            [0,     1,      -4.5,   4.5,    0,      4.99e-16,      1.2767e-15,      0 ,     -4.996e-16,      0],
                            [0,     0,      0,      0,      9.,     -22.5,  13.5,   -22.5,  27.,    13.5],
                            [0,     0,      0,      0,      0,      27.,    -27.,   0,      -27.,   0],
                            [0,     0,      0,      0,      0,      -4.5,   13.5,   0,      0,      0],
                            [0,     0,      0,      0,      -4.5,   4.5,    0 ,     18.,    -13.5,  -13.5],
                            [0,     0,      0,      0,      0,      -4.5,   0,      0,      13.5,   0],
                            [0,     0,      0,      0,      1,      4.7184e-16,      -6.3837e-16,      -4.5,   -1.2767e-15,      4.5]
                        ])
    j = 0
    for i in range(10):
        j=0
        for s in range(3):
            for r in range(3):
                phi[i] += coef[i,j]*(xi**r)*(eta**s)
                if(r!=0):
                    dphi[i,0] += coef[i,j]*(xi**(r-1))*(eta**s)*r
                if(s!=0):
                    dphi[i,1] += coef[i,j]*(xi**r)*(eta**(s-1))*s
                j = j+1        
    return phi,dphi


def quad1d(q):
    # Order 5 Gauss-Legendre points
    if(q == 3 ):
        x_v = np.array([
            0.112701665379258, 0.500000000000000, 0.887298334620742
        ])
        w_v = np.array([
            0.277777777777778, 0.444444444444444, 0.277777777777778
        ])
    # Order 7 Gauss-Legendre points
    elif(q == 4 ):
        x_v = np.array([
            0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026
        ])
        w_v = np.array([
            0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
        ])
    #  Order 9 Gauss-Legendre points
    elif(q == 5 ):
        x_v = np.array([
            0.046910077030668, 0.230765344947158, 0.500000000000000, 0.769234655052841, 0.953089922969332
        ])
        w_v = np.array([
            0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095
        ])
     #  Order 11 Gauss-Legendre points
    elif(q == 6 ):
        x_v = np.array([
            0.033765242898424, 0.169395306766868, 0.380690406958402, 0.619309593041598,
            0.830604693233132, 0.966234757101576
        ])
        w_v = np.array([
            0.085662246189585, 0.180380786524069, 0.233956967286345, 0.233956967286345,
            0.180380786524069, 0.085662246189585
        ])
    return x_v, w_v

def quad2d():
    x_v = np.array([
            0.333333333333333, 0.333333333333333, 0.020634961602525, 0.489682519198738,
            0.489682519198738, 0.489682519198738, 0.489682519198738, 0.020634961602525,
            0.125820817014127, 0.437089591492937, 0.437089591492937, 0.437089591492937,
            0.437089591492937, 0.125820817014127, 0.623592928761935, 0.188203535619033,
            0.188203535619033, 0.188203535619033, 0.188203535619033, 0.623592928761935,
            0.910540973211095, 0.044729513394453, 0.044729513394453, 0.044729513394453,
            0.044729513394453, 0.910540973211095, 0.036838412054736, 0.221962989160766,
            0.221962989160766, 0.741198598784498, 0.741198598784498, 0.036838412054736,
            0.221962989160766, 0.036838412054736, 0.741198598784498, 0.221962989160766,
            0.036838412054736, 0.7411985987844
        ])
    w_v = np.array([
            0.048567898141400, 0.015667350113570, 0.015667350113570, 0.015667350113570,
            0.038913770502387, 0.038913770502387, 0.038913770502387, 0.039823869463605,
            0.039823869463605, 0.039823869463605, 0.012788837829349, 0.012788837829349,
            0.012788837829349, 0.021641769688645, 0.021641769688645, 0.021641769688645,
            0.021641769688645, 0.021641769688645, 0.021641769688645
        ])
    x_v = x_v.reshape((19,2))
    return x_v, w_v



def face_node(V, local_face,s):
    if(local_face==0):
        V = np.delete(V,[0,1,2,4,5,7],0)
    elif(local_face==1):
        V = np.delete(V,[1,2,3,5,6,8],0)
        V = V[::-1,:]
    elif(local_face==2):
        V = np.delete(V,[4,5,6,7,8,9],0)
    phi,_ = basis_geom_1d(s)
    x = np.sum(phi*V[:,0])
    y = np.sum(phi*V[:,1])
    return x,y

def Jacobian_edge(V, local_face,s):
    if(local_face==0):
        V = np.delete(V,[0,1,2,4,5,7],0)
    elif(local_face==1):
        V = np.delete(V,[1,2,3,5,6,8],0)
        V = V[::-1,:]
    elif(local_face==2):
        V = np.delete(V,[4,5,6,7,8,9],0)
    phi,dphi = basis_geom_1d(s)
    J1 = np.sum(dphi*V[:,0])
    J2 = np.sum(dphi*V[:,1])
    return J1,J2

def Jacobian_nonlinear(V, xi, eta):
    _, dphi = basis_geom_2d(xi, eta)
    J1 = np.sum(dphi[:,0]*V[:,0])
    J2 = np.sum(dphi[:,1]*V[:,0])
    J3 = np.sum(dphi[:,0]*V[:,1])
    J4 = np.sum(dphi[:,1]*V[:,1])
    det  = J1*J4 - J2*J3 
    return det

def Area_curved(V, npp):
    M = np.zeros((npp,npp))
    sq,wq = quad2d()
    for q in range(19):
        xi = sq[q,0]; eta = sq[q,1]
        phi = basis_2d(npp,xi,eta)
        det = Jacobian_nonlinear(V,xi, eta)
        for i in range(npp):
            for j in range(npp):
                M[i,j] += phi[i]*phi[j]*det*wq[q]
    A = np.sum(M)
    return A 