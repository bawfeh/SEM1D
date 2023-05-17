# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
from matplotlib import pyplot as plt

inputfile = 'LinAdv0002'

filename = './InputFile/'+inputfile+'.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

# # ============================================
# 
# import pickle
# with open('BurgersGRKDG4p6.pickle', 'rb') as f:
# # with open('BurgersGRKDG4p6t1.pickle', 'rb') as f:
#     exactsol = pickle.load(f)
# from timeStepping import ExpSL 
# expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
# nstp = len(exactsol[0])-1
# expsl.uc = exactsol[0][int(nstp/4)]
# # ===========================================

import Initdata
bc = Initdata.Source(inputfile)

from timeStepping import SLRKDG as Method
# from timeStepping import RKDG as Method

N = 6
err = np.zeros((N,))
dt = np.zeros((N,))
PP =  np.zeros((N,),dtype=int)
Cr = np.zeros((N,))
sol = None
# U = 0.75
U = 1.
# U = 2.75
# nsteps = 80

for n in range(N):

    gll = GaussMesh.GaussDistribution(data.polyDeg)
    mesh = GaussMesh.MeshDG(data, gll)
    
#     data.Cr = 1./mesh.getNe()
    lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
    nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
    midpoints = (nodes[1:]+nodes[:-1])/2
    dx = min(midpoints[1:]-midpoints[:-1])
#     dx = mesh.size(0)
    ds = data.Cr*dx / U
    nsteps = max(int(data.dT/ds),1)
    print("Cr = %f, Ne = %i." %(data.Cr, mesh.getNe()))
#     ds = data.dT / nsteps
#     Cr[n] = U * ds / dx
#     print("Cr = %f, Ne = %i." %(Cr[n], mesh.getNe()))

    timestepping = Method(data, gll, mesh, bc, nsteps)
#     timestepping.stepBystepReport = False 

    timestepping.evolve()

    if data.flowtype['nonlinear']:
        expsl.nv = mesh.nodes.size
        sol = expsl.InterpSL(mesh.nodes) # exact sol on numerical grid
    else: sol = timestepping.init.compute4ex(data.tb[-1])

    err[n] = timestepping.errorL2(timestepping.sol[-1], 0, sol) if data.storeAll\
            else timestepping.errorL2(timestepping.sol, 0, sol)
    dt[n] = timestepping.dt
    PP[n] = data.polyDeg
    print("After %i steps (q = %i), L2-Error = %8.4E." %(nsteps,PP[n],err[n]))
#     plt.plot(mesh.nodes, timestepping.sol,'--', mesh.nodes, sol)
#     plt.show()
    data.polyDeg += 1
#     data.Cr /= 2

    if data.flowtype['nonlinear']:
        normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
        err[n] /= normsol  # relative error
# ===========================================

plt.semilogy(PP, err,'.-')
plt.xlabel('q')
plt.ylabel('Error (L2)')
plt.show()

if data.storeAll:
    plt.plot(mesh.nodes,timestepping.sol[-1],mesh.nodes,sol,'--')
else:
    plt.plot(mesh.nodes,timestepping.sol,mesh.nodes,sol,'--')
plt.xlabel('x')
plt.show()

# ============================================
with open('temp1c.npy', 'wb') as f:
    np.save(f, err)
    np.save(f, PP)
    np.save(f, dt)
# ===========================================

slopes = (np.log(err[1:])-np.log(err[:-1])) / (PP[1:]-PP[:-1])
print('Errors =',err)
print('slopes =',slopes)
print('Ne =', mesh.getNe())


# with open('temp1a.npy','rb') as f:
#     err = np.load(f)
#     PP = np.load(f)
#     dt = np.load(f)
# errs = np.zeros((3,len(err)))
# dT = np.zeros((3,len(dt)))
# errs[0] = err
# dT[0] = dt
# with open('temp1b.npy','rb') as f:
#     err = np.load(f)
# errs[1] = err
# dT[1] = dt
# with open('temp1c.npy','rb') as f:
#     err = np.load(f)
# errs[2] = err
# dT[2] = dt
# with open('LinearRKDG1pConvCr2p5.npy', 'wb') as f:
#     np.save(f, errs)
#     np.save(f, PP)
#     np.save(f, dT)
# 
# plt.semilogy(PP, errs[0], PP, errs[1], PP, errs[2])
# plt.xlabel('q')
# plt.ylabel(r'Error (L$_2$)')
# plt.show()
