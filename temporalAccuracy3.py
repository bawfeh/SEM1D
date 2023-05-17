# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh

filename = '/Users/bawfeh78/Documents/Cplusplus/SLDG/input5.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GaussMesh.GaussDistribution(data.polyDeg)
mesh = GaussMesh.MeshDG(data, gll)
    
U = .75
# U = 1.
lnext = mesh.nodes[1:].tolist()+[mesh.nodes[0]]
nodes = mesh.nodes[abs(np.array(lnext)-mesh.nodes) > 1e-15]
midpoints = (nodes[1:]+nodes[:-1])/2
dx = min(midpoints[1:]-midpoints[:-1])
#     dx = mesh.size(0)
ds = data.Cr*dx / U
nsteps = int(data.dT/ds)
# nsteps = 125
print('Cr =', U * (data.dT/nsteps) / dx)

# ============================================

import pickle
# with open('BurgersGRKDG3.pickle', 'rb') as f:
# with open('BurgersGRKDG4p6t1.pickle', 'rb') as f:
with open('BurgersGRKDG4p6.pickle', 'rb') as f:
    exactsol = pickle.load(f)
from timeStepping import ExpSL 
expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
nstp = len(exactsol[0])-1
expsl.uc = exactsol[0][int(nstp/2)]
expsl.nv = mesh.nodes.size
sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid
# ===========================================

import Initdata as Initd
bc = Initd.Source('NlinAdv0005d',data.Re)

from timeStepping import SLRKDG
timestepping = SLRKDG(data, gll, mesh, bc, nsteps)
# from timeStepping import RKDG
# timestepping = RKDG(data, gll, mesh, bc, nsteps)
# timestepping.stepBystepReport = False 


N = 10
err = np.zeros((N,))
dt = np.zeros((N,))
Cr =  np.zeros((N,))

for n in range(N):

    timestepping.evolve()
    err[n] = timestepping.errorL2(timestepping.sol[-1], 0, sol) if data.storeAll\
            else timestepping.errorL2(timestepping.sol, 0, sol)
    dt[n] = timestepping.dt
    Cr[n] = U * dt[n] / dx
    print("After %i steps (Cr = %f) L2-Error = %f." %(timestepping.nsteps, Cr[n], err[n]))
    data.Cr += 0.5
    ds = data.Cr*dx / U
    nsteps = int(data.dT/ds)
    timestepping.nsteps = nsteps
#     timestepping.nsteps *= 2
    timestepping.dt = data.dT / timestepping.nsteps
    if n > 0:
        order = (np.log(err[n])-np.log(err[n-1])) / (np.log(dt[n])-np.log(dt[n-1]))
        print('order =', order)


normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
err /= normsol  # relative error
# ===========================================

# with open('temp4.npy','wb') as f:
#     np.save(f, dt)
#     np.save(f, err)
#     np.save(f, Cr)

from matplotlib import pyplot as plt

# if data.storeAll:
#     plt.plot(mesh.nodes,timestepping.sol[-1],mesh.nodes,sol)
# else:
#     plt.plot(mesh.nodes,timestepping.sol,mesh.nodes,sol)
# plt.xlabel('x')
# plt.show()

h1, = plt.loglog(dt, err,'.-')
h2, = plt.loglog(dt, dt**data.orderRK,'--')
plt.xlabel('h')
plt.ylabel('Error (L2)')
plt.legend([h1,h2],['SLRKDG%i'%data.orderRK,'h^%i'%data.orderRK])
plt.show()

plt.semilogy(Cr, err,'.-')
plt.xlabel(r'$\sigma$')
plt.ylabel(r'Error (L$_2$)')
plt.show()
# # ============================================
slopes = (np.log(err[1:])-np.log(err[:-1])) / (np.log(dt[1:])-np.log(dt[:-1]))

print('slopes =',slopes)
print('Cr =', Cr)
# 
# # ===========================================
# errs = []
# CCr = []
with open('temp1a.npy','wb') as f:
    np.save(f, err)
    np.save(f, dt)
    np.save(f, Cr)
# errs.append(err)
# CCr.append(Cr)
# with open('temp3.npy','rb') as f:
#     dt = np.load(f)
#     err = np.load(f)
#     Cr = np.load(f)
# errs.append(err)
# CCr.append(Cr)
# with open('BurgersGRKDG3xB3.npy','rb') as f:
#     dt = np.load(f)
#     err = np.load(f)
#     Cr = np.load(f)
# errs.append(err)
# CCr.append(Cr)
# with open('BurgersGRKDG3temp.npy','wb') as ff:
#     np.save(ff, errs)
#     np.save(ff, dt)
#     np.save(ff, CCr)
# ===========================================
