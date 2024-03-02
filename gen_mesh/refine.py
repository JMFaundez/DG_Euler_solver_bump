"""
Code to refine mesh
"""

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

# read files of mesh that will be refined.
file_num = 0
print("Refining mesh", file_num)
I2E = np.loadtxt('I2E'+str(file_num)+'.txt') - 1 #python index
B2E = np.loadtxt('B2E'+str(file_num)+'.txt') - 1 #python index
mesh = readgri('bump'+str(file_num)+'.gri')
nodes = mesh['V']
elements = mesh['E']
boundaries = mesh['B']
N_elem = len(elements)

df_nodes = pd.DataFrame({'x':nodes[:,0],'y':nodes[:,1]})
df_I2E = pd.DataFrame({'elemL':I2E[:,0],'faceL':I2E[:,1],'elemR':I2E[:,2],'faceR':I2E[:,3], 'node1':"", 'node2':"","new_node":""})
df_B2E = pd.DataFrame({'elem':B2E[:,0],'face':B2E[:,1],'bgroup':B2E[:,2], 'node1':"", 'node2':"","new_node":""})
N_face = len(I2E)
N_nodes = len(nodes)
N_boundaries = len(B2E)

index_new_nodes = N_nodes - 1 #index number to store in file
i_p = [1,2,0] #mapping for the next local face. Ej: face[ip[0]]=1
i_m = [2,0,1] #mapping for the previous local face. Ej: face[im[0]]=2

#loop over the interior faces
for i in range(N_face):
    elem = int(df_I2E.loc[i,"elemL"])
    index = int(df_I2E.loc[i,"faceL"])
    n1 = int(elements[elem][i_p[index]])
    n2 = int(elements[elem][i_m[index]])
    #store nodes of the face
    df_I2E.at[i,'node1'] = n1
    df_I2E.at[i,'node2'] = n2
    index_new_nodes = index_new_nodes + 1 # number of the new node
    #store number of the new node
    df_I2E.at[i,'new_node'] = index_new_nodes
    #get coordinates of the nodes of the face
    x1 = df_nodes.loc[n1,"x"]
    y1 = df_nodes.loc[n1,"y"]
    x2 = df_nodes.loc[n2,"x"]
    y2 = df_nodes.loc[n2,"y"]
    #nx and ny are the coordinates for the new node
    nx = (x1+x2)/2
    ny = (y1+y2)/2
    #store coordinates of the new node (the number is the index of the dataframe)
    df_newnode = pd.DataFrame(data={"x":[nx],"y":[ny]})
    #df_nodes = df_nodes.append(df_newnode,ignore_index=True)
    df_nodes = pd.concat([df_nodes,df_newnode],ignore_index=True)

#loop over the boundary faces
for i in range(N_boundaries):
    elem = int(df_B2E.loc[i,"elem"])
    index = int(df_B2E.loc[i,"face"])
    bgroup = int(df_B2E.loc[i,"bgroup"])
    n1 = int(elements[elem][i_p[index]])
    n2 = int(elements[elem][i_m[index]])
    df_B2E.at[i,'node1'] = n1
    df_B2E.at[i,'node2'] = n2
    x1 = df_nodes.loc[n1,"x"]
    y1 = df_nodes.loc[n1,"y"]
    x2 = df_nodes.loc[n2,"x"]
    y2 = df_nodes.loc[n2,"y"]
    nx = (x1+x2)/2
    #print('x2:',x2)
    #print('nx:',nx)
    #check if the vboudnary is bottom
    if bgroup == 0:
        ny = 0.0625*np.exp(-25*nx**2)
        df_nodes.loc[n1,"y"] = 0.0625*np.exp(-25*x1**2)
        df_nodes.loc[n2,"y"] = 0.0625*np.exp(-25*x2**2)
    else:
        ny = (y1+y2)/2
    #store coordinates of the new node
    df_newnode = pd.DataFrame(data={"x":[nx],"y":[ny]})
    #df_nodes = df_nodes.append(df_newnode,ignore_index=True)
    df_nodes = pd.concat([df_nodes,df_newnode],ignore_index=True)
    index_new_nodes = index_new_nodes + 1
    df_B2E.at[i,'new_node'] = index_new_nodes

#plt.figure()
#plt.plot(df_nodes['x'],df_nodes['y'],'.')
#plt.axis('equal')
#plt.show()
#list of the new elements
new_elements = []
for i in range(N_elem):
    node_e = elements[i] #old nodes of the element
    face_l = df_I2E.loc[df_I2E.elemL==i] #left faces of the element
    face_r = df_I2E.loc[df_I2E.elemR==i] #right faces of the element
    face_b = df_B2E.loc[df_B2E.elem==i] #boundary faces of the element
    num_l = len(face_l.index) #number of left faces
    num_r = len(face_r.index) #number of right faces
    num_b = len(face_b.index) #number of boudnary faces
    center_element = set([]) #this is the element created by the new nodes, the center triangle
    #for each node is created a new element
    for j in range(3):
        n_new = [] #new element for each node
        nodej = node_e[j]
        #each node has to be in two faces, if it's found ==> store the new node of that face
        #the new element is created using the old node + the new nodes of the others faces
        for k in range(num_l):
            if face_l.iloc[k]["node1"] == nodej or face_l.iloc[k]["node2"] == nodej:
                n_new.append(face_l.iloc[k]["new_node"])
                center_element.add(face_l.iloc[k]["new_node"])

        for k in range(num_r):
            if face_r.iloc[k]["node1"] == nodej or face_r.iloc[k]["node2"] == nodej:
                n_new.append(face_r.iloc[k]["new_node"])
                center_element.add(face_r.iloc[k]["new_node"])

        for k in range(num_b):
            if face_b.iloc[k]["node1"] == nodej or face_b.iloc[k]["node2"] == nodej:
                n_new.append(face_b.iloc[k]["new_node"])
                center_element.add(face_b.iloc[k]["new_node"])
        #append the old node to the new element
        n_new.append(nodej)
        #append the new element associated to the old node
        new_elements.append(n_new)
    #append the 'center' element
    center_element = list(center_element)
    new_elements.append(center_element)

new_elements = np.array(new_elements)+1 #1-index file
#store elements and nodes, change name according to the mesh
fileelem = 'element'+str(file_num+1)
np.savetxt(fileelem,new_elements,fmt='%i')
print("File ",fileelem," saved")
new_nodes = df_nodes.values
filenodes='nodes'+str(file_num+1)
np.savetxt(filenodes,new_nodes)
print("File ",filenodes," saved")