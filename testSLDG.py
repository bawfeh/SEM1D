# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
import Initdata
import time

inputfile = 'LinAdv0001b'

filename = './InputFile/'+inputfile+'.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GaussMesh.GaussDistribution(data.polyDeg)

mesh = GaussMesh.MeshDG(data, gll)

bc = Initdata.Source(inputfile)
init = Initdata.Initdata(bc, mesh)

u0 = init.compute4ex(data.tb[-1])

# mesh.grid()
# mesh.plot(u0)

# nsteps = 256
# U = 2.75
# U = .75
U = 1.
# lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
# nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
# midpoints = (nodes[1:]+nodes[:-1])/2
# dx = min(midpoints[1:]-midpoints[:-1])
# dt = data.Cr*dx / U
dt = data.Cr * data.dx[0] / U
# for constant velocity, we choose Cr integer, and dt = Cr*h 
# this permits xtics arrive at grid points
nsteps = int(data.dT/dt)
print('h =',data.dx[0])

sol = []
normsol = 1
nsub = 5
# for k in range(nsteps, 1, -1):
#     if k%(nsub-1) == 0: nsteps = k; break

# nsteps = 15360
# nsteps = 8042
# nsteps = 80
dt = data.dT/nsteps
print('dx = %5.2f' % mesh.size(0))
print('Cr = %5.2f; q = %i; nsteps = %i; dt = %3.2fh = %9.5f' % (U * dt/data.dx[0], data.polyDeg, nsteps, dt/data.dx[0],dt))


# from timeStepping import SLRKDG
# timestepping = SLRKDG(data, gll, mesh, bc, nsteps)
from timeStepping import RKDG
timestepping = RKDG(data, gll, mesh, bc, nsteps)
# timestepping.exactsol = True
# timestepping.stepBystepReport = False 

tic = time.time()

timestepping.evolve()

toc = time.time()

print("Exection time:", toc-tic)


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
#             normsol = timestepping.errorL2(np.zeros(np.shape(sol1)), 0, sol1) 
            err[k] = timestepping.errorL2(timestepping.sol[ind[k]], 0, sol1)# /normsol 
    else:
        for k in range(nsub):
            sol1 = timestepping.init.compute4ex(t[k])
            sol.append(sol1)
            normsol = timestepping.errorL2(np.zeros(np.shape(sol1)), 0, sol1) 
            err[k] = timestepping.errorL2(timestepping.sol[ind[k]], 0, sol1)# /normsol 
    print("L2-errors at times", t, ": err = ", err)
else:
    t = timestepping.nstepsdone * timestepping.dt
    if data.flowtype['nonlinear']: 
        expsl.uc = exactsol[0][-1]
        sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
    else:
        sol = timestepping.init.compute4ex(t)
#     normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
    err = timestepping.errorL2(timestepping.sol, 0, sol) 
#     err /= normsol  # relative error
    print("L2-Error obtain after %i step(s) = %8.4E." %(timestepping.nstepsdone,err / timestepping.nstepsdone)) # test test test

# ===========================================

from matplotlib import pyplot as plt

if data.storeAll:
    ind =  np.linspace(0, nsteps, nsub, dtype=int)
    for k in range(nsub):
        plt.plot(mesh.nodes,timestepping.sol[ind[k]])
        plt.plot(mesh.nodes,sol[k],'--')
else:
    plt.plot(mesh.nodes,sol)
    plt.plot(mesh.nodes,timestepping.sol,'r--')
    for el in mesh.elements:
        plt.plot(mesh.nodes[el],timestepping.sol[el],'g.-')
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
