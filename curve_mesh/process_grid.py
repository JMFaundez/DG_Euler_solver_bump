"""
 Code to read .gri file and compute 
 matrices for task2. It also calculates
 the normal for task3. 
"""

#%%
import numpy as np
import pandas as pd


def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)


#-----------------------------------------------------------
#read .gri (1-based index) and convert to 0-based index
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

#returns the local face number (0,1 or 2)
def node_face(n1,n2,element):
    element = np.asarray(element)
    res = np.isin(element,[n1,n2])
    res = [not i for i in res]
    return int(np.where(res)[0])

#compute the normal vector given two points using anticlockwise numbering
def normal_vec(n1,n2,nodes):
    x1,y1 = nodes[int(n1)]
    x2,y2 = nodes[int(n2)]
    dx = x2-x1 
    dy = y2-y1
    l = np.sqrt((x2-x1)**2+(y2-y1)**2)
    return dy/l, -dx/l


#compute area for a given element
def compute_area(element,nodes):
    p0 = nodes[int(element[0])]
    p1 = nodes[int(element[1])]
    p2 = nodes[int(element[2])]
    suma = (p0[0]*p1[1] - p0[1]*p1[0]) + (p1[0]*p2[1] - p1[1]*p2[0]) + (p2[0]*p0[1] - p2[1]*p0[0])  
    return np.abs(suma)/2

# Generate the matrices for the mesh
def get_matrices(file_name):
    #read file
    mesh = readgri(file_name)
    #nx and ny are the x and y componentes of the normal vector
    #node1 and node2 are the numbers of the nodes of the face
    d_I2E = {'elemL':[], 'faceL':[], 'elemR':[], 'faceR':[], 'node1':[], 'node2':[], 'nx':[],'ny':[]}
    dfI2E = pd.DataFrame(data=d_I2E)
    d_B2E = {'elem':[], 'face':[], 'bgroup':[], 'nx':[],'ny':[]}
    dfB2E = pd.DataFrame(data=d_B2E)
    nodes = mesh['V']
    boundaries = mesh['B']
    Ec = mesh['Ec']
    El = mesh['El']
    columns = [1, 2, 4, 5 , 6, 7, 8]
    Ec = np.delete(Ec, columns, 1)
    elements = np.concatenate((El,Ec))
    N_elem = len(elements)
    Area = np.zeros((N_elem,1))

    #loop over each element of the mesh
    for i in range(N_elem):
        n1, n2, n3 = elements[i]
        Area[i] = compute_area(elements[i],nodes)
        faces  = [[n1,n2],[n2,n3],[n3,n1]]
        for j in range(3):
            facej = faces[j] #nodes/faces of the element
            # check if facej is a boundary face
            is_boundary = False #flag to break loop
            for bo in range(4): #4 boundaries
                boundary = boundaries[bo]
                for ib in range(len(boundary)): #loop over each boundary
                    #boolean values to check if face belong to boundaries
                    #nodes can be clockwise or anticlockwise
                    con1 = (facej[0] == boundary[ib][0] and facej[1] == boundary[ib][1])
                    con2 = (facej[1] == boundary[ib][0] and facej[0] == boundary[ib][1]) 
                    if con1 or con2:
                        is_boundary = True #change flag
                        #compute normal vector
                        nx, ny = normal_vec(facej[0],facej[1],nodes)
                        #store face and normal vector to dataframe
                        dfb = pd.DataFrame({'elem':[i], 'face':[node_face(facej[0],facej[1],elements[i])], 'bgroup':[bo], 'nx':[nx],'ny':[ny]})
                        #dfB2E = dfB2E.append(dfb, ignore_index='True')
                        dfB2E = pd.concat([dfB2E,dfb])
                        break
                if is_boundary:
                    break
            
            if is_boundary:
                continue
            #if face isn't boundary:
            else:
                #check if the face was evaluated, if so it is the right element
                if len(dfI2E.loc[(dfI2E.node1==facej[1]) & (dfI2E.node2==facej[0])])==1:
                    dfI2E.loc[(dfI2E.node1==facej[1]) & (dfI2E.node2==facej[0]),'elemR'] = i
                    dfI2E.loc[(dfI2E.node1==facej[1]) & (dfI2E.node2==facej[0]),'faceR'] = node_face(facej[0],facej[1],elements[i])
                #check if the left face is empty
                elif len(dfI2E.loc[(dfI2E.node1==facej[0]) & (dfI2E.node2==facej[1])])==0:
                    nx, ny = normal_vec(facej[0],facej[1],nodes)
                    df1 = pd.DataFrame(data={'elemL':[i], 'faceL':[node_face(facej[0],facej[1],elements[i])], 'elemR':[0], 'faceR':[0], 'node1':[facej[0]], 'node2':[facej[1]],'nx':[nx],'ny':[ny]})
                    #dfI2E = dfI2E.append(df1, ignore_index='True')
                    dfI2E = pd.concat([dfI2E,df1])
    interior = dfI2E.values #convert the dataframe to array
    I2E = interior[:,0:4]+1 #1-index based matrix
    In = interior[:,6:]
    I2E_nodes = interior[:,4:6]
    bound = dfB2E.values
    B2E = bound[:,0:3]+1 #1-index based matrix
    Bn = bound[:,3:]
    result = {'I2E_nodes':I2E_nodes,'I2E':I2E,'B2E':B2E,'In':In,'Bn':Bn,'Area':Area}
    return result


#read file
file_num = '1'
file_name = 'bump'+file_num+'_curv.gri'
mesh = readgri(file_name)
nodes = mesh['V']
boundaries = mesh['B']
Ec = mesh['Ec']
El = mesh['El']
columns = [1, 2, 4, 5 , 6, 7, 8]
Ec = np.delete(Ec, columns, 1)
elements = np.concatenate((El,Ec))
print('# nodes: ',nodes.shape,'# elements', elements.shape)
N_elem = len(elements)
sum_normal = np.zeros(N_elem) #initilaze vector of the sum of normal vector
#actually do task2
T2 = get_matrices(file_name)
#store I2E and B2E matrices, change name according to mesh
np.savetxt('I2E'+file_num+'.txt',T2['I2E'],fmt='%i')
np.savetxt('B2E'+file_num+'.txt',T2['B2E'],fmt='%i')
np.savetxt('In'+file_num+'.txt',T2['In'],fmt='%.14f')
np.savetxt('Bn'+file_num+'.txt',T2['Bn'],fmt='%.14f')
np.savetxt('Area'+file_num+'.txt',T2['Area'],fmt='%.14f')
print('# Boundary ',T2['Bn'].shape,'# Interior ',T2['In'].shape)
I2E = T2['I2E'] - 1 # 0-index
B2E = T2['B2E'] - 1 #0-index
N_in = len(I2E) #number of interior faces
N_bn = len(B2E) #number of boundary faces
In = np.zeros((N_in,2)) # normal vector interior faces
Bn = np.zeros((N_bn,2)) #normal vector boundary faces
Area = T2['Area']
i_p = [1,2,0] #mapping for the next local face. Ej: face[ip[0]]=1
i_m = [2,0,1] #mapping for the  previous local face. Ej: face[im[0]]=2
sum_n = np.zeros((N_elem,2))
mag_el = np.zeros(N_elem)

#compute (normal_vec)&length for each interior face
for i in range(N_in):
    ele_1 = elements[int(I2E[i,0])]
    ele_2 = elements[int(I2E[i,2])]
    n1 = ele_1[i_p[int(I2E[i,1])]]
    n2 = ele_2[i_p[int(I2E[i,3])]]
    nx,ny = normal_vec(n1,n2,nodes)
    dx = nodes[int(n2)][0] - nodes[int(n1)][0]
    dy = nodes[int(n2)][1] - nodes[int(n1)][1]
    l = np.sqrt(dx**2+dy**2)
    In[i,0] = nx*l
    In[i,1] = ny*l

#compute (normal_vec)&length for each boundary face
for i in range(N_bn):
    ele = elements[int(B2E[i,0])]
    n1 = ele[i_p[int(B2E[i,1])]]
    n2 = ele[i_m[int(B2E[i,1])]]
    nx,ny = normal_vec(n1,n2,nodes)
    dx = nodes[int(n2)][0] - nodes[int(n1)][0]
    dy = nodes[int(n2)][1] - nodes[int(n1)][1]
    l = np.sqrt(dx**2+dy**2)
    Bn[i,0] = nx*l
    Bn[i,1] = ny*l

#get the faces for each element and sum outward normal vectors
for i in range(N_elem):
    index_b = np.isin(B2E[:,0], i) #find the boundary faces that belong to the i element
    index_il = np.isin(I2E[:,0], i) #find the left interior faces that belong to the i element
    index_ir = np.isin(I2E[:,2], i) #find the right nterior faces that belong to the i element
    bound_vec = Bn[index_b,:]
    left_vec = In[index_il,:]
    right_vec = In[index_ir,:]
    su_b = np.sum(bound_vec,0)
    su_l = np.sum(left_vec,0)
    su_r = -np.sum(right_vec,0)
    sum_n[i,:] = su_b+su_l+su_r
    mag_el[i] = np.sqrt(sum_n[i,0]**2+sum_n[i,1]**2)

print(max(mag_el))








            
        


    
  

            
