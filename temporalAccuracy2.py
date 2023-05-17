# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh

inputfile = 'LinAdv0005'

filename = './InputFile/'+inputfile+'.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GaussMesh.GaussDistribution(data.polyDeg)

mesh = GaussMesh.MeshDG(data, gll)

# nsteps = 256
# U = .75
# U = 2.
U = 1.
# U = 10.
lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
midpoints = (nodes[1:]+nodes[:-1])/2
dx = min(midpoints[1:]-midpoints[:-1])
dt = data.Cr*dx / U
nsteps = int(data.dT/dt)

# nsteps = 15360
# nsteps = 8042
# ============================================

import Initdata
bc = Initdata.Source(inputfile)
init = Initdata.Initdata(bc, mesh)

import pickle
# with open('BurgersGRKDG4p6.pickle', 'rb') as f:
# with open('NCF2RKDG3.pickle', 'rb') as f:
with open('BLRKDG3.pickle', 'rb') as f:
    exactsol = pickle.load(f)
from timeStepping import ExpSL 
expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
expsl.nv = mesh.nodes.size
# nsub = 21
# nsub = 9
nsub = 5
t = None
nstp = len(exactsol[0])-1
sol = np.zeros((nsub,expsl.nv)) if data.storeAll else np.zeros((expsl.nv,))
if data.storeAll:
    ind = np.linspace(0, nstp, nsub, dtype=int)
    t = np.linspace(data.tb[0],data.tb[-1],nstp+1)[ind]
    print('exact solutions at times t =', t)
    if data.flowtype['nonlinear']: 
        for k in range(nsub):
            expsl.uc = exactsol[0][ind[k]]
            sol[k] = expsl.InterpSL(mesh.nodes) # exact sol on numerical grid
    else:
        for k in range(nsub):
            sol[k] = init.compute4ex(t[k])
else:
    if data.flowtype['nonlinear']: 
        expsl.uc = exactsol[0][-1]
        sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
    else: sol = init.compute4ex(data.tb[-1])
# ===========================================

# for k in range(nsteps, nstp):
for k in range(nsteps, 1, -1):
    if nstp%k == 0 and k%(nsub-1) == 0: nsteps = k; break

dt = data.dT/nsteps
print('Cr =', U * dt/dx, '; nsteps =', nsteps, '; dt =', dt)

from timeStepping import SLRKDG
timestepping = SLRKDG(data, gll, mesh, bc, nsteps)
# from timeStepping import RKDG
# timestepping = RKDG(data, gll, mesh, bc, nsteps)
# timestepping.stepBystepReport = False 

timestepping.evolve()

err = np.zeros((nsub,)) if data.storeAll else 0.
ind =  np.linspace(0, nsteps, nsub, dtype=int)
if data.storeAll:
    t = np.linspace(data.tb[0],data.tb[-1],nsteps+1)[ind]
    for k in range(nsub):
        normsol = timestepping.errorL2(np.zeros(np.shape(sol[k])), t[k], sol[k]) 
        err[k] = timestepping.errorL2(timestepping.sol[ind[k]],t[k], sol[k]) / normsol
    print("L2-Errors obtained at times", t, " =", err)
else:
    normsol = timestepping.errorL2(np.zeros(np.shape(sol)), data.tb[-1], sol)
    err = timestepping.errorL2(timestepping.sol, data.tb[-1], sol) / normsol  
    print("L2-Error obtain at time t = %f steps = %5.4E." %(data.tb[-1], err))
# # ===========================================
# 
# with open('BLRKDG3xB.npy', 'wb') as f:
#     np.save(f, np.array(timestepping.sol)[ind])
#     np.save(f, mesh.nodes)
#     np.save(f, t)
#     np.save(f, sol)
#     np.save(f, err)
# # ===========================================

from matplotlib import pyplot as plt

if data.storeAll:
    for k in range(nsub):
        plt.plot(mesh.nodes,timestepping.sol[ind[k]])
        plt.plot(mesh.nodes,sol[k],'--')
else:
    plt.plot(mesh.nodes,sol)
    plt.plot(mesh.nodes,timestepping.sol,'r.--')
    plt.xlabel('x')

# umin = -0.25
# umax = 0.75
# plt.plot(mesh.nodes, umax*np.ones(np.shape(mesh.nodes)), 'k--')
# plt.plot(mesh.nodes, umin*np.ones(np.shape(mesh.nodes)), 'k--')
plt.show()

# if data.storeAll:
#     plt.semilogy(t[1:], err[1:], marker='.')
#     plt.show()
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
