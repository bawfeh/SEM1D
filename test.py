# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
import Initdata
import time

inputfile = 'NLinAdv0005c'

filename = './InputFile/'+inputfile+'.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GaussMesh.GaussDistribution(data.polyDeg)

mesh = GaussMesh.Mesh(data, gll)

bc = Initdata.Source(inputfile, data.Re)
init = Initdata.Initdata(bc, mesh)

u0 = init.compute4ex()

# mesh.grid()
# mesh.plot(u0)


# nsteps = 64
nsteps = 64
U = 1.
lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
midpoints = (nodes[1:]+nodes[:-1])/2
dx = min(midpoints[1:]-midpoints[:-1])
# # dt = data.Cr*dx / U
# dt = data.Cr * data.dx[0]
# # for constant velocity, we choose Cr integer, and dt = Cr*h 
# # this permits xtics arrive at grid points
# nsteps = int(data.dT/dt)
print('h =',data.dx[0])
dt = data.dT/nsteps
# print('Cr =', U * dt/dx, '; nsteps =', nsteps, '; dt =', dt)
print("Cr = %4.2f; nsteps = %i; dt = %f = %4.1fdx" %(U * dt/dx, nsteps, dt, dt/data.dx[0]))

from timeStepping import CF_Asch
timestepping = CF_Asch(data, gll, mesh, bc, nsteps)
# from timeStepping import IMEX_Asch
# timestepping = IMEX_Asch(data, gll, mesh, bc, nsteps)

start = time.time()
timestepping.evolve2()
print("Execution time:", time.time()-start)

# mesh.plot(timestepping.sol)

# ============================================
import pickle
with open('odeRadau1.pickle', 'rb') as f:
    exactsol = pickle.load(f)

from timeStepping import ExpSL 

expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
expsl.uc = np.zeros((expsl.nv,))
expsl.nv = mesh.nodes.size
l = 0
for el in expsl.mesh.elements: # scatter operation (copy exactsol[0])
    expsl.uc[expsl.mesh.elementsL[l]] = exactsol[0][:,-1][el]; l += 1
sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
# ===========================================
from matplotlib import pyplot as plt

if data.storeAll:
    ntot = len(timestepping.sol)
    for k in range(ntot): 
        if np.isnan(timestepping.sol[k]).any(): continue
        if np.isinf(timestepping.sol[k]).any(): continue
        plt.plot(mesh.nodes,timestepping.sol[k])
        plt.show()
    print('L2-error =', timestepping.errorL2(timestepping.sol[-1], 0, sol))
else:
    h1, = plt.plot(mesh.nodes,sol,'r')
    h2, = plt.plot(mesh.nodes,timestepping.sol,'.')
    plt.xlabel('x')
    plt.legend([h1,h2],['exact','SLRK3'])
    print('L2-error =', timestepping.errorL2(timestepping.sol, 0, sol))
    plt.show()
# ============================================
