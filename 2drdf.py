#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#default parameters
npoints = 100  #numbers of shells
length = 15   #length maximum of pair distance 

#read coordination from POSCAR
f = open("POSCAR", "r+")
tmp = f.readlines()
atom_num = int(tmp[6])
axis = []
for i in [2,3,4]:
    axis.append(list(map(float, tmp[i].split())))
axis = np.array(axis)
scale_factor = float(tmp[1])

atoms = []
for i in range(8,8+atom_num):
    atoms.append(list(map(float, tmp[i].split())))
atoms = np.array(atoms)
f.close()

#get area of cell
tmp = np.array([axis[0,1]*axis[1,2]-axis[0,2]*axis[1,1],axis[0,2]*axis[1,0]-\
        axis[0,0]*axis[1,2],axis[0,0]*axis[1,1]-axis[0,1]*axis[1,0]])
area = np.linalg.norm(tmp)
area *= scale_factor**2

#calculating rdf
g = np.zeros(npoints)
delta = length/npoints #distance between adjent shells

for i in range(atom_num):
    for j in range(i+1,atom_num):
        a = atoms[i].tolist()
        b = atoms[j].tolist()
        for k in range(3):
            if a[k] - b[k]> 0.5:
                a[k] -= 1
            if a[k] - b[k]<-0.5:
                b[k] -= 1
        a = np.dot(a,axis)*scale_factor
        b = np.dot(b,axis)*scale_factor
        tmp = a-b
        tmp = np.linalg.norm(tmp[:2])
        num = int(tmp/delta)
        if num < npoints:
            g[num] += 2

#averaging
rho = atom_num/area
for i in range(npoints):
    g[i] /= 2*np.pi*(i+1)*delta*delta*rho

#plot
x = [delta*(i+1) for i in range(npoints)]
plt.plot(x,g)
plt.show()
