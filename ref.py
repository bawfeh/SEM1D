# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
import Initdata
import time

inputfile = 'NLinAdv0005d'
outputf = 'BurgersGRKDG'

filename = './InputFile/'+inputfile+'.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GaussMesh.GaussDistribution(data.polyDeg)

mesh = GaussMesh.MeshDG(data, gll)

bc = Initdata.Source(inputfile)
init = Initdata.Initdata(bc, mesh)

u0 = init.compute4ex(data.tb[0])

# mesh.grid()
# mesh.plot(u0)

# nsteps = 256
# U = 10.
# U = 2.75
U = .75
# U = 1.

lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
midpoints = (nodes[1:]+nodes[:-1])/2

dx = min(midpoints[1:]-midpoints[:-1])
# dx = data.dx[0]

nsteps = 15360
# nsteps = 8042
# nsteps = 320

nsub = 5

if data.storeAll:
    if np.mod(nsteps, (nsub-1)) != 0:
        nsb = nsub-1
        nsteps = nsb * ((nsteps // nsb) + 1)

dt = data.dT/nsteps

print('Cr = %5.2f; q = %i; nsteps = %i; dt = %3.2fh = %6.5f; T = %2.2f'
        % (U * dt/dx, data.polyDeg, nsteps, dt/data.dx[0],dt, data.dT))

# from timeStepping import SLRKDG
# timestepping = SLRKDG(data, gll, mesh, bc, nsteps)
from timeStepping import RKDG
timestepping = RKDG(data, gll, mesh, bc, nsteps)
# timestepping.exactsol = True
# timestepping.stepBystepReport = False 

# ===========================================
tic = time.time()

timestepping.evolve()

print("Exection time:", time.time()-tic)
# ===========================================


# import pickle
# with open(outputf+'4p6.pickle', 'wb') as f:
#     pickle.dump([timestepping.sol, gll, mesh], f)
# ===========================================


print('Cr = %5.2f; q = %i; nsteps = %i; Ne = M = %i; dt = %3.2fh = %6.5f; T = %2.2f'
        % (U * dt/dx, data.polyDeg, nsteps, data.Ne, dt/data.dx[0],dt, data.dT))

from matplotlib import pyplot as plt

if data.storeAll:
    ind =  np.linspace(0, nsteps, nsub, dtype=int)
    for k in range(nsub):
        plt.plot(mesh.nodes,timestepping.sol[ind[k]])
    plt.xlabel('x')
else:
    plt.plot(mesh.nodes,timestepping.sol,'r--')
    plt.xlabel('x')

plt.show()
# ===========================================
