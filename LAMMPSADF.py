#!/usr/bin/python3
#this script create file0 and file1 for LAMMPS ADF calculation
#ponychen
#20190722
#email:18709821294@outlook.com

import numpy as np

orientation1 = list(map(float, input("write in your first crystal orientation:  ").split()))
orientation2 = list(map(float, input("write in your second crystal orientation:  ").split()))
orientation3 = list(map(float, input("write in your third crystal orientation: ").split()))
lattice=float(input("please input your fcc axis constant:  "))

orientation = np.array([orientation1, orientation2, orientation3])

for i in range(3):
    length = sum(orientation[i,:]**2)**0.5
    for j in range(3):
        orientation[i,j] /= length

atoms = np.array([[0.5,0.5,0],[-0.5,0.5,0],[0.5,0,-0.5],[-0.5,0,-0.5],\
        [0,0.5,-0.5],[0,-0.5,-0.5]])

coor = []
for i in range(6):
    a = lattice*np.dot(orientation, atoms[i,:])
    coor.append(a.tolist())

f = open("fcc.vec", "a+")
for data in coor:
    line = map(str, data)
    line = "    ".join(line)
    line += "\n"
    f.write(line)
f.close()
