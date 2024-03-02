import numpy as np

def new_coord(c, x1,y1,x2,y2):
    l = np.sqrt((x1-x2)**2+(y1-y2)**2)
    xn = x1 + c*(x2-x1)
    yn = y1 + c*(y2-y1)
    return xn, yn

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

# read files of mesh that will be curved.
q = 3; # it's the only value support for now
file_num = 1
B2E = np.loadtxt('../gen_mesh/B2E'+str(file_num)+'.txt') - 1 #python index
mesh = readgri('../gen_mesh/bump'+str(file_num)+'.gri')
nodes = mesh['V']
elements = mesh['E']
boundaries = mesh['B']
names = mesh['Bname']
N_elem = len(elements)
N_bound = len(B2E)

index_boundaries = []
local_faces = []
for i in range(N_bound):
    if(B2E[i,2]==0):
        index_boundaries.append(int(B2E[i,0]))
        local_faces.append(int(B2E[i,1]))

el_bot = elements[index_boundaries]
N_bot = len(el_bot)

nq = int((q+1)*(q+2)/2)
new_elem = np.zeros((N_bot,nq))

N_new_nodes = (nq-3)*N_bot
new_nodes = np.zeros((N_new_nodes,2))
nod = 0
if(q==3):
    for i in range(N_bot):
        elem = index_boundaries[i]
        face = local_faces[i]
        n1 = elements[elem,0]
        n2 = elements[elem,1]
        n3 = elements[elem,2]
        new_elem[i,0] = n1
        new_elem[i,3] = n2
        new_elem[i,9] = n3
        x1 = nodes[n1,0]
        y1 = nodes[n1,1]
        x2 = nodes[n2,0]
        y2 = nodes[n2,1]
        x3 = nodes[n3,0]
        y3 = nodes[n3,1]
        if(face==0):
            # Node 1
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(1./3, x1,y1,x2,y2)
            new_elem[i,1] = nod + len(nodes)
            nod +=1 
            # Node 2
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(2./3, x1,y1,x2,y2)
            new_elem[i,2] = nod + len(nodes)
            nod +=1
            #Node 3
            nodes[n2,1] =  0.0625*np.exp(-25.*x2**2)
            # Node 4
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(2./3, x3,y3,x1,y1)
            new_elem[i,4] = nod + len(nodes)
            nod +=1
            # Node 5
            new_nodes[nod,0] = (x1+x2+x3)/3
            new_nodes[nod,1] = (y1+y2+y3)/3
            new_elem[i,5] = nod + len(nodes)
            nod +=1 
            # Node 6 //python index
            xn, _ = new_coord(1./3,x2,y2,x3,y3)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,6] = nod + len(nodes)
            nod +=1 
            # Node 7
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(1./3, x3,y3,x1,y1)
            new_elem[i,7] = nod + len(nodes)
            nod +=1
            # Node 8
            xn,_ = new_coord(2./3,x2,y2,x3,y3)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,8] = nod + len(nodes)
            nod +=1
            #Node 9
            nodes[n3,1] =  0.0625*np.exp(-25.*x3**2)
        elif(face==1):
            # Node 0
            nodes[n1,1] =  0.0625*np.exp(-25.*x1**2)
            # Node 1
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(1./3, x1,y1,x2,y2)
            new_elem[i,1] = nod + len(nodes)
            nod +=1 
            # Node 2
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(2./3, x1,y1,x2,y2)
            new_elem[i,2] = nod + len(nodes)
            nod +=1
            # Node 4
            xn,_ = new_coord(2./3, x3,y3,x1,y1)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,4] = nod + len(nodes)
            nod +=1
            # Node 5
            new_nodes[nod,0] = (x1+x2+x3)/3
            new_nodes[nod,1] = (y1+y2+y3)/3
            new_elem[i,5] = nod + len(nodes)
            nod +=1 
            # Node 6 //python index
            new_nodes[nod,0], new_nodes[nod,1]= new_coord(1./3,x2,y2,x3,y3)
            new_elem[i,6] = nod + len(nodes)
            nod +=1 
            # Node 7
            xn, _ = new_coord(1./3, x3,y3,x1,y1)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,7] = nod + len(nodes)
            nod +=1
            # Node 8
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(2./3,x2,y2,x3,y3)
            new_elem[i,8] = nod + len(nodes)
            nod +=1
            # Node 9
            nodes[n3,1] =  0.0625*np.exp(-25.*x3**2)
        elif(face==2):
            # Node 0
            nodes[n1,1] =  0.0625*np.exp(-25.*x1**2)
            # Node 1
            xn, _ = new_coord(1./3, x1,y1,x2,y2)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,1] = nod + len(nodes)
            nod +=1 
            # Node 2
            xn, _ = new_coord(2./3, x1,y1,x2,y2)
            yn = 0.0625*np.exp(-25.*xn**2)
            new_nodes[nod,0] = xn 
            new_nodes[nod,1] = yn
            new_elem[i,2] = nod + len(nodes)
            nod +=1
            # Node 0
            nodes[n2,1] =  0.0625*np.exp(-25.*x2**2)
            # Node 4
            new_nodes[nod,0], new_nodes[nod,1]=  new_coord(2./3, x3,y3,x1,y1)
            new_elem[i,4] = nod + len(nodes)
            nod +=1
            # Node 5
            new_nodes[nod,0] = (x1+x2+x3)/3
            new_nodes[nod,1] = (y1+y2+y3)/3
            new_elem[i,5] = nod + len(nodes)
            nod +=1 
            # Node 6 //python index
            new_nodes[nod,0], new_nodes[nod,1]= new_coord(1./3,x2,y2,x3,y3)
            new_elem[i,6] = nod + len(nodes)
            nod +=1 
            # Node 7
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(1./3, x3,y3,x1,y1)
            new_elem[i,7] = nod + len(nodes)
            nod +=1
            # Node 8
            new_nodes[nod,0], new_nodes[nod,1] = new_coord(2./3,x2,y2,x3,y3)
            new_elem[i,8] = nod + len(nodes)
            nod +=1
            
f = open('bump'+str(file_num)+'_curv.gri',"w+")
nNode = len(nodes) + len(new_nodes)
Dim = 2 
nElemTot = len(elements)
#write first line
f.write("%i %i %i\n" % (nNode,nElemTot, Dim))
#write old nodes
for i in range(len(nodes)):
    f.write("%.16f %.16f\n" % (nodes[i,0],nodes[i,1]))
#write new nodes
for i in range(len(new_nodes)):
    f.write("%.16f %.16f\n" % (new_nodes[i,0],new_nodes[i,1]))

nBGroup = 4
#write number of boundaries
f.write("%i\n" % (nBGroup))

for i in range(nBGroup):
    nBFace = len(boundaries[i])
    f.write("%i %i " % (nBFace,2)+names[i]+"\n")
    for j in range(nBFace):
        f.write("%i %i\n" % (boundaries[i][j][0]+1,boundaries[i][j][1]+1))

# write linear elements
elements_lin = np.delete(elements,index_boundaries,0)
f.write("%i %i " % (len(elements_lin),1)+"TriLagrange\n")
for i in range(len(elements_lin)):
    f.write("%i %i %i\n" % (elements_lin[i][0]+1,elements_lin[i][1]+1,elements_lin[i][2]+1))

# write curved elements
new_elem = new_elem +1
f.write("%i %i " % (len(new_elem),q)+"TriLagrange\n")
for i in range(len(new_elem)):
    f.write("%i %i %i %i %i %i %i %i %i %i\n" % (new_elem[i][0], new_elem[i][1], new_elem[i][2], new_elem[i][3], new_elem[i][4], new_elem[i][5], new_elem[i][6], new_elem[i][7], new_elem[i][8], new_elem[i][9]))