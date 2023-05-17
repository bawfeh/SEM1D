
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
# ============================================
import pickle
with open('odeRadau1.pickle', 'rb') as f:
    exactsol = pickle.load(f)

from timeStepping import ExpSL 
expsl = ExpSL(exactsol[-1],exactsol[1],data.xb)
expsl.uc = np.zeros((expsl.nv,))
expsl.nv = mesh.nodes.size
# ===========================================

from timeStepping import CF_Asch
# from timeStepping import IMEX_Asch
import Initdata as Initd

method = 'SL233'

Cr = 1.8
midpoints = (mesh.nodes[1:]+mesh.nodes[0:-1])/2
dx = min(midpoints[1:]-midpoints[0:-1])
print('dx',dx)
dt = Cr*dx
nsteps = int(data.dT/dt)
N = 10
nu = 1e-3 * np.linspace(1,N,N)
# timestepping.stepBystepReport = False 
err = np.zeros((N,))

for n in range(N):
    l = 0
    for el in expsl.mesh.elements: # scatter operation (copy exactsol[0])
        expsl.uc[expsl.mesh.elementsL[l]] = exactsol[0][:,n][el]; l += 1
    sol = expsl.InterpSL(mesh.nodes) # exact sol or numerical grid

    data.Re = nu[n]
    bc = Initdata.Source(inputfile, data.Re)
    timestepping = CF_Asch(data, gll, mesh, bc, nsteps)
#     timestepping = IMEX_Asch(data, gll, mesh, bc, nsteps)
    timestepping.evolve()
    err[n] = timestepping.errorL2(timestepping.sol, 0, sol) 
    print("L2-Error obtain with nu = %f is given by error = %f." %(nu[n],err[n]))

    normsol = timestepping.errorL2(np.zeros(np.shape(sol)), 0, sol) 
    err[n] /= normsol  # relative error
# ============================================

print(err)

# with open('varyViscosity'+method+'.npy','wb') as f:
#     np.save(f, err)
#     np.save(f, nu)
#     np.save(f, dt)
#     np.save(f, nsteps)
# # ===========================================
# from matplotlib import pyplot as plt
# 
# h1, = plt.plot(mesh.nodes,timestepping.sol)
# h2, = plt.plot(mesh.nodes,sol)
# plt.xlabel('x')
# plt.legend([h1,h2],['num','exact'])
# plt.show()
# 
# plt.loglog(nu, err,'.-')
# plt.xlabel('nu')
# plt.ylabel('Error (L2)')
# plt.show()
# # ============================================
# # 
# # with open('temp1.npy','rb') as f:
# #     err = np.load(f)
# # errs = np.zeros((len(err),3))
# # errs[:,0] = err
# # with open('temp2.npy','rb') as f:
# #     err = np.load(f)
# # errs[:,1] = err
# # with open('temp3.npy','rb') as f:
# #     err = np.load(f)
# #     nu = np.load(f)
# # errs[:,2] = err
# # with open('varyViscosityIMEX.npy','wb') as f:
# #     np.save(f, errs)
# #     np.save(f, nu)
# # 
