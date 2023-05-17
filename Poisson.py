# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile as IN
import GaussMesh as GM

filename = '/Users/bawfeh78/Documents/Cplusplus/SLDG/input1.txt'

data = IN.Parameters() # default input data

flag = IN.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GM.GaussDistribution(data.polyDeg)

mesh = GM.Mesh(data, gll)

import Assembly
from LinearSolvers import LinearSolvers

semAss =  Assembly.Assembly(mesh,gll)

lsolve = LinearSolvers(semAss,1e-8,2000)

# ===========================================

x = mesh.nodes

uex = np.cos(np.pi * x)

data.Dirbc_data = uex[mesh.Dnodes]

semAss.setbc(data.Dirbc_data)

fs = np.pi**2 * np.cos(np.pi * x)

# lsolve.setSolverType('biCGstab',True)
# lsolve.setSolverType('biCGstab')
# lsolve.setSolverType('CG')

if (lsolve.matvec):
    f = semAss.scatter(fs)
    semAss.Mass(f)
    f = semAss.dssum(f)
else: f = fs

u = lsolve.solve(f)

# ===========================================
from matplotlib import pyplot as plt

if (lsolve.matvec): u = semAss.gather(u, False, False)

plt.plot(x,u,x,uex,'r--')

plt.show()
