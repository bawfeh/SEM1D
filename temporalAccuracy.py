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
method = 'SL233'

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

nsteps = 10
from timeStepping import CF_Asch
timestepping = CF_Asch(data, gll, mesh, bc, nsteps)
# from timeStepping import IMEX_Asch
# timestepping = IMEX_Asch(data, gll, mesh, bc, nsteps)
# timestepping.stepBystepReport = False 
N = 8
err = np.zeros((N,))
dt = np.zeros((N,))

start = time.time()
for n in range(N):
    timestepping.evolve()
    err[n] = timestepping.errorL2(timestepping.sol, 0, sol) 
    dt[n] = timestepping.dt
    print("L2-Error obtain with %i steps = %f." %(nsteps,err[n]))
    nsteps *= 2
    timestepping.nsteps = nsteps
    timestepping.nsteps2do = nsteps
    timestepping.dt = data.dT / nsteps

print("Execution time:", time.time()-start)
normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
err /= normsol  # relative error
# ===========================================
from matplotlib import pyplot as plt

plt.plot(mesh.nodes,timestepping.sol,mesh.nodes,sol)
plt.xlabel('x')
plt.show()

h1, = plt.loglog(dt, err,'.-')
h2, = plt.loglog(dt, dt**3,'--')
plt.xlabel('h')
plt.ylabel('Error (L2)')
plt.legend([h1,h2],['SLRK3G','h^3'])
plt.show()
# # ============================================
# with open('temp1.npy', 'rb') as f:
#     errs = np.load(f)
# # errs = np.zeros((len(err),3))
# errs[:,1] = err
# with open('temp2.npy','wb') as f:
#     np.save(f, errs)
#     np.save(f, dt)
# ===========================================
slopes = (np.log(err[1:])-np.log(err[:-1])) / (np.log(dt[1:])-np.log(dt[:-1]))

print('errors =',err)
print('slopes =',slopes)

# with open('temporalAccuracy'+method+'.npy','wb') as f:
#     np.save(f, err)
#     np.save(f, dt)
#     np.save(f, data.Ne)
#     np.save(f, data.polyDeg)
#     np.save(f, data.Re)
#     np.save(f, data.tb)
# 
