# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
from matplotlib import pyplot as plt

inputfile = 'LinAdv0004'

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
# expsl.uc = exactsol[0][-1]#[int(nstp/2)]
# # ===========================================

import Initdata as Initd
bc = Initd.Source(inputfile, data.Re)

from timeStepping import SLRKDG
from timeStepping import RKDG
# 
N = 6
err = np.zeros((N,))
dt = np.zeros((N,))
Ne =  np.zeros((N,),dtype=int)
sol = None
# U = 0.75
# U = 2.75
U = 1.
L = data.xb[-1] - data.xb[0]

for n in range(N):

    gll = GaussMesh.GaussDistribution(data.polyDeg)
    mesh = GaussMesh.MeshDG(data, gll)
    
#     lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
#     nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
#     midpoints = (nodes[1:]+nodes[:-1])/2
#     dx = min(midpoints[1:]-midpoints[:-1])
    dx = mesh.size(0)
    ds = data.Cr*dx / U
    nsteps = max(10,int(data.dT/ds))
    print("Cr = %f, Ne = %i, dt = %3.1fh = %6.4f." %(U*ds/dx, mesh.getNe(), ds*data.Ne/L, ds))

    timestepping = SLRKDG(data, gll, mesh, bc, nsteps)
#     timestepping = RKDG(data, gll, mesh, bc, nsteps)
#     timestepping.stepBystepReport = False 

    if data.flowtype['nonlinear']: 
        expsl.nv = mesh.nodes.size
        sol = expsl.InterpSL(mesh.nodes) # exact sol on numerical grid
    else: sol = timestepping.init.compute4ex(data.tb[-1])

    timestepping.evolve()
    err[n] = timestepping.errorL2(timestepping.sol[-1], 0, sol) if data.storeAll\
            else timestepping.errorL2(timestepping.sol, 0, sol)
    dt[n] = timestepping.dt
    Ne[n] = mesh.getNe()
    print("After %i steps (Ne = %i), L2-Error = %5.4E." %(nsteps,Ne[n],err[n]))# / timestepping.nstepsdone ))
#     plt.plot(mesh.nodes, sol, mesh.nodes, timestepping.sol,'--')
#     plt.show()
    data.Ne  *= 2
    data.nx[0] = data.nx[0] * 2
    data.dx[0] = data.dx[0] / 2

#     data.Cr /= 2
# 
#     normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
#     err[n] /= normsol  # relative error
#     err[n] /= timestepping.nstepsdone 
    if (n > 0):
        slope = np.log(err[n]/err[n-1]) / np.log(Ne[n-1]/Ne[n])
        print("order =", slope)
# ===========================================
# from matplotlib import pyplot as plt


# h1, = plt.loglog(dt, err,'.-')
# h2, = plt.loglog(dt, dt**data.orderRK,'--')
# plt.xlabel('h')
# plt.ylabel('Error (L2)')
# plt.legend([h1,h2],['SLRK3G','h^%i'%data.orderRK])
# plt.show()

# L = data.xb[-1]-data.xb[0]

plt.loglog(L/Ne, err,'.-')
plt.xlabel('h')
plt.ylabel(r'Error (L$_2$)')
plt.show()

if data.storeAll:
    plt.plot(mesh.nodes,timestepping.sol[-1],mesh.nodes,sol,'--')
else:
    plt.plot(mesh.nodes,timestepping.sol,mesh.nodes,sol,'--')
plt.xlabel('x')
plt.show()

# # ============================================
# with open('temp1a.npy', 'wb') as f:
#     np.save(f, err)
#     np.save(f, Ne)
#     np.save(f, dt)
# # ===========================================

print('Ne =', Ne)
print('L =', L)

slopes = (np.log(err[1:])-np.log(err[:-1])) / (np.log(Ne[:-1])-np.log(Ne[1:]))

print('errors =',err)
print('slopes =',slopes)
# 
# with open('temp1a.npy','rb') as f:
#     err = np.load(f)
#     Ne = np.load(f)
#     dt = np.load(f)
# errs = np.zeros((5,len(err)))
# errs[0] = err
# 
# with open('temp1b.npy','rb') as f:
#     err = np.load(f)
# errs[1] = err
# 
# with open('temp1c.npy','rb') as f:
#     err = np.load(f)
# errs[2] = err
# 
# with open('temp1d.npy','rb') as f:
#     err = np.load(f)
# errs[3] = err
# 
# with open('temp1e.npy','rb') as f:
#     err = np.load(f)
# errs[4] = err
# 
# with open('LinearRKDG3Cr3p5Convp3.npy', 'wb') as f:
#     np.save(f, errs)
#     np.save(f, Ne)
#     np.save(f, dt)
# 
# L = 2*np.pi #data.xb[-1]-data.xb[0]
# with open('LinearRKDG3Cr3p5Convp3.npy', 'rb') as f:
#     errs = np.load(f)
#     Ne = np.load(f)
# plt.loglog(L/Ne, errs[0], L/Ne, errs[1], L/Ne, errs[2], L/Ne, errs[3], L/Ne, errs[4] )
# plt.xlabel('2/M')
# plt.ylabel(r'Error (L$_2$)')
# plt.show()
# print('Ne =', Ne)
# 
# slopes = (np.log(errs[:,1:])-np.log(errs[:,:-1])) / (np.log(Ne[:-1])-np.log(Ne[1:]))
# print('slopes =',slopes)
