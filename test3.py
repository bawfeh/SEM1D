# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
import Initdata
import time

# from timeStepping import SLRKDG as Method
from timeStepping import RKDG as Method


inputfile = 'NLinAdv0006b'
outputf = 'NCFRKDG'

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
# U = 2.75
# U = .75
U = 10.
# U = 1.

lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
midpoints = (nodes[1:]+nodes[:-1])/2

dx = min(midpoints[1:]-midpoints[:-1])
# dx = data.dx[0]

# for constant velocity, we choose Cr integer, and dt = Cr*h 
# this permits xtics depart from grid points

# dt = data.Cr*dx / U
# nsteps = max(2,int(data.dT/dt))
# nsteps = 80
nsteps = 160
print('h =',data.dx[0])

sol = []
normsol = 1
nsub = 5

if data.storeAll:
    if np.mod(nsteps, (nsub-1)) != 0:
        nsb = nsub-1
        nsteps = nsb * ((nsteps // nsb) + 1)

dt = data.dT/nsteps
print('Cr = %5.2f; q = %i; nsteps = %i; dt = %3.2fh = %6.5f; T = %2.2f'
        % (U * dt/dx, data.polyDeg, nsteps, dt/data.dx[0],dt, data.dT))

# ============================================
import pickle
with open('./'+outputf+'/'+outputf+'4p6.pickle', 'rb') as f:
    exactsol = pickle.load(f)
from timeStepping import ExpSL 
expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
expsl.nv = mesh.nodes.size
nstp = len(exactsol[0])-1
lnd = np.linspace(0, nstp, nsub, dtype=int)
# ===========================================

timestepping = Method(data, gll, mesh, bc, nsteps)
# timestepping.exactsol = True
# timestepping.stepBystepReport = False 

tic = time.time()

timestepping.evolve()

print("Exection time:", time.time()-tic)
# ===========================================

err = None
if data.storeAll and (timestepping.nstepsdone==timestepping.nsteps) :
    ind =  np.linspace(0, nsteps, nsub, dtype=int)
    err = np.zeros((nsub,))
    t = np.linspace(data.tb[0], data.tb[-1], nsub)
    if data.flowtype['nonlinear']: 
        for k in range(nsub):
            expsl.uc = exactsol[0][lnd[k]]
            sol1 = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
            sol.append(sol1)
            normsol = timestepping.errorL2(np.zeros(np.shape(sol1)), 0, sol1) 
            err[k] = timestepping.errorL2(timestepping.sol[ind[k]], 0, sol1) /normsol 
    else:
        for k in range(nsub):
            sol1 = timestepping.init.compute4ex(t[k])
            sol.append(sol1)
#             normsol = timestepping.errorL2(np.zeros(np.shape(sol1)), 0, sol1) 
            err[k] = timestepping.errorL2(timestepping.sol[ind[k]], 0, sol1)# /normsol 
    print("Rel.L2-errors at times", t, ": err = ", err)
else:
    t = timestepping.nstepsdone * timestepping.dt
    if data.flowtype['nonlinear']: 
        expsl.uc = exactsol[0][-1]
        sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
        err = timestepping.errorL2(timestepping.sol, 0, sol) 
        normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
        err /= normsol  # relative error
    else:
        sol = timestepping.init.compute4ex(t)
        err = timestepping.errorL2(timestepping.sol, 0, sol) 
#         normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
#         err /= normsol  # relative error
    print("Rel. L2-Error obtain after %i step(s) = %8.4E." %(timestepping.nstepsdone,err))

# ===========================================
# 
# ind =  np.linspace(0, nsteps, nsub, dtype=int)
# nsol = []; t = np.linspace(data.tb[0], data.tb[-1], nsub)
# for k in ind:
#     nsol.append(timestepping.sol[k])
# 
# with open(outputf+'3SL.npy', 'wb') as f:
#     np.save(f, nsol)
#     np.save(f, mesh.nodes)
#     np.save(f, t)
#     np.save(f, sol)
#     np.save(f, err)
# ===========================================

print('Cr = %5.2f; q = %i; nsteps = %i; Ne = M = %i; dt = %3.2fh = %6.5f; T = %2.2f'
        % (U * dt/dx, data.polyDeg, nsteps, data.Ne, dt/data.dx[0],dt, data.dT))

from matplotlib import pyplot as plt

if data.storeAll:
    ind =  np.linspace(0, nsteps, nsub, dtype=int)
    for k in range(nsub):
        plt.plot(mesh.nodes,timestepping.sol[ind[k]])
        plt.plot(mesh.nodes,sol[k],'--')
else:
    plt.plot(mesh.nodes,sol)
    plt.plot(mesh.nodes,timestepping.sol,'r--')
#     for el in mesh.elements:
#         plt.plot(mesh.nodes[el],timestepping.sol[el],'g.-')
    plt.xlabel('x')

plt.show()
# ===========================================

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
# 
# if data.storeAll:
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     x = mesh.nodes
#     ntot = len(timestepping.sol)
#     t = np.arange(0,ntot) * dt + data.tb[0]
#     X, Y = np.meshgrid(t, x)
#     Z = np.array(timestepping.sol).T
#     surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
#     fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()
# # ============================================
