
#%%
import matplotlib.pyplot as plt
import numpy as np 
import time

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
    # read elements
    Ne0 = 0; E = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0==0) else np.concatenate((E,Ei), axis=0)
        Ne0 += ne
    f.close()
    Mesh = {'V':V, 'E':E, 'B':B, 'Bname':Bname }
    return Mesh

filen="1"
file_name = "bump"+filen+".gri"
faces = np.loadtxt("I2E"+filen+".txt")-1
bfaces = np.loadtxt("B2E"+filen+".txt")-1
mesh = readgri(file_name)
p = mesh['V']
t = mesh['E']
linew=1
N_bfaces = len(bfaces[:,0])
N_faces = len(faces[:,0])

t_i = time.time()
#plt.axis('off')
plt.figure()
plt.ylim([0,1])
i_p = [1,2,0] #mapping for the next local face. Ej: face[ip[0]]=1
i_m = [2,0,1]
for i in range(N_faces):
    e1,ix,e2,_ = faces[i]
    n1 = int(t[int(e1),i_p[int(ix)]])
    n2 = int(t[int(e1),i_m[int(ix)]])
    x1 = p[n1,0]
    y1 = p[n1,1]
    x2 = p[n2,0]
    y2 = p[n2,1]
    plt.axis('equal')
    plt.plot([x1,x2],[y1,y2],'k-',lw=linew)

l_s = ['r-','g-','b-','y-']
for i in range(N_bfaces):
    e1,ix,nb = bfaces[i]
    n1 = int(t[int(e1),i_p[int(ix)]])
    n2 = int(t[int(e1),i_m[int(ix)]])
    x1 = p[n1,0]
    y1 = p[n1,1]
    x2 = p[n2,0]
    y2 = p[n2,1]
    plt.axis('equal')
    plt.plot([x1,x2],[y1,y2],l_s[int(nb)],lw=linew)
plt.tight_layout
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("mesh"+filen+".png", bbox_inches='tight',pad_inches=0)
t_f = time.time()
delta_t = (t_f-t_i)/60
