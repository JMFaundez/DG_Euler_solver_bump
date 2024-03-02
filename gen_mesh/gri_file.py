"""
 Code to generate .gri file 
"""
#%%
import numpy as np 
import matplotlib.pyplot as plt
filetosave="bump0.gri"
f = open(filetosave,"w+")

p = np.loadtxt('nodes')
t = np.loadtxt('element')
tol = 0.02 # Tolerance to find boundaries
nNode, Dim = np.shape(p)
nElemTot,_ = np.shape(t)
print("Number of elements:",nElemTot)
print("Number of nodes:",nNode)
#write first line
f.write("%i %i %i\n" % (nNode,nElemTot, Dim))

#write node coordinates
for i in range(nNode):
        if np.abs(p[i,1] - 0.0625*np.exp(-25.*p[i,0]**2))<tol:
                f.write("%.16f %.16f\n" % (p[i,0],0.0625*np.exp(-25.*p[i,0]**2)))
        else:
                f.write("%.16f %.16f\n" % (p[i,0],p[i,1]))

#Number of boundaries
nBGroup = 4
#write number of boundaries
f.write("%i\n" % (nBGroup))


index_left = np.where((np.abs(p[:,0] +1.5))<tol)
p_left = p[index_left]
index_left = np.asarray(index_left)
index_left_sort = index_left[0,np.argsort(p_left[:,1])]

index_right = np.where(np.abs(p[:,0] - 1.5)<tol)
p_right = p[index_right]
index_right = np.asarray(index_right)
index_right_sort = index_right[0,np.argsort(p_right[:,1])]

index_top = np.where(np.abs(p[:,1] -0.8)<tol)
p_top = p[index_top]
index_top = np.asarray(index_top)
index_top_sort = index_top[0,np.argsort(p_top[:,0])]

index_bot = []
for i in range(nNode): #check bottom nodes
    if np.abs(p[i,1] - 0.0625*np.exp(-25.*p[i,0]**2))<tol:
        index_bot.append(i)
        
p_bot = p[index_bot]
index_bot = np.asarray(index_bot)
index_bot_sort = index_bot[np.argsort(p_bot[:,0])]
print("========== Number of nodes at each boundary =========== ")
print("#Nodes left:",len(index_top[0,:]))
print("#Nodes right:",len(index_right[0,:]))
print("#Nodes left:",len(index_left[0,:]))
print("#Nodes bottom:",len(index_bot))
Group_names = ['Bottom', 'Right', 'Top', 'Left']
#nodes are sorted to save them
index_sort = [index_bot_sort,index_right_sort, index_top_sort, index_left_sort]

for i in range(nBGroup):
    nBFace = len(index_sort[i])-1
    f.write("%i %i " % (nBFace,2)+Group_names[i]+"\n")
    for j in range(nBFace):
        f.write("%i %i\n" % (index_sort[i][j]+1,index_sort[i][j+1]+1))

f.write("%i %i " % (nElemTot,1)+"TriLagrange\n")

#function to check if the nodes of the elements are sorted clockwise/counter
def clockwise(element):
    p0 = p[int(element[0]-1)]
    p1 = p[int(element[1]-1)]
    p2 = p[int(element[2]-1)]
    suma = (p0[0]*p1[1] - p0[1]*p1[0]) + (p1[0]*p2[1] - p1[1]*p2[0]) + (p2[0]*p0[1] - p2[1]*p0[0])  
    return suma

for i in range(nElemTot):
    if clockwise(t[i])>0: #check if the nodes are saved counterclockwise
        f.write("%i %i %i\n" % (t[i][0],t[i][1],t[i][2]))
    else:
        f.write("%i %i %i\n" % (t[i][0],t[i][2],t[i][1]))

f.close()
print("File ",filetosave," saved!")


