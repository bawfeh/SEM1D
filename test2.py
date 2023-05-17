# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
import Initdata
import time

inputfile = 'NLinAdv0007'

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

nsteps = 2**8
midpoints = (mesh.nodes[1:]+mesh.nodes[:-1])/2
u0 = init.compute4Ex(mesh.nodes[:-1], data.tb[0])
Ds = midpoints[1:]-midpoints[:-1]
Ds = np.append(Ds,midpoints[0]-mesh.nodes[0]+mesh.nodes[-1]-midpoints[-1])
dx = min(Ds)
print('dx =',dx)
# dt = Cr / max(abs(u0/Ds))
# print('dt =',dt)
# nsteps = int(data.dT/dt)
# print('nsteps =', nsteps)
dt = data.dT/nsteps
Cr = dt * max(abs(u0/Ds))
# print('Cr =', Cr)
print("Cr = %4.2f; nsteps = %i; dt = %f = %10.6fdx" %(Cr, nsteps, dt, dt/data.dx[0]))

from timeStepping import CF_Asch
timestepping = CF_Asch(data, gll, mesh, bc, nsteps)
# from timeStepping import IMEX_Asch
# timestepping = IMEX_Asch(data, gll, mesh, bc, nsteps)
timestepping.evolve()

# from timeStepping import ExactSol
# 
# from scipy.integrate import solve_ivp
# 
# exactsol = ExactSol(data, gll, mesh, bc)
# sol = solve_ivp(exactsol.feval, [exactsol.t0, exactsol.tend], exactsol.u0, method='Radau',atol=1e-8, rtol=1e-8)
# 
# import pickle
# with open('odeRadauKdV05c.pickle', 'wb') as f:
#     pickle.dump([sol.y[:,-1], gll, mesh], f)
# 
# ===========================================
from matplotlib import pyplot as plt
# import pickle
# with open('odeRadauKdV.pickle', 'rb') as f:
# # with open('odeRadauKdV05.pickle', 'rb') as f:
#     exactsol = pickle.load(f)

# # err = np.linalg.norm(timestepping.sol - exactsol[0], np.inf) / np.linalg.norm(exactsol[0], np.inf)
# # print('errinf =', err)
# # sol = exactsol[0][:,-1]
# ln = len(exactsol[0][:,-1])
# # sol = exactsol[0][np.arange(0,ln,10)]
# ref = timestepping.errorLinf(0., 0, exactsol[0][:,-1])

if data.storeAll:
    ntot = len(timestepping.sol)
    for k in range(ntot): 
        if np.isnan(timestepping.sol[k]).any(): continue
        if np.isinf(timestepping.sol[k]).any(): continue
        plt.plot(mesh.nodes,timestepping.sol[k])
        plt.show()
#     print('L2-error =', timestepping.errorLinf(timestepping.sol[-1], 0, sol) / ref)
else:
    h1, = plt.plot(mesh.nodes,timestepping.sol)
#     h2, = plt.plot(exactsol[-1].nodes, exactsol[0][:,-1])
#     plt.plot(mesh.nodes,sol.y[:,-1])
    plt.xlabel('x')
#     plt.legend([h1,h2],['num','ref'])
#     print('L2-error =', timestepping.errorLinf(timestepping.sol, 0, sol) / ref)
    plt.show()
# ============================================
