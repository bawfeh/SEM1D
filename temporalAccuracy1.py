# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile as IN
import GaussMesh as GM
from scipy.optimize import least_squares

filename = '/Users/bawfeh78/Documents/Cplusplus/SLDG/input2.txt'

data = IN.Parameters() # default input data

flag = IN.InputFile(filename).updateData(data) # with file inputs

if flag : print('\nData successfully updated!\n')
else: print('\nWARNING: Incorrect data provided in '+ filename)

gll = GM.GaussDistribution(data.polyDeg)

mesh = GM.Mesh(data, gll)

# ============================================
import pickle
with open('odeRadauKdV05.pickle', 'rb') as f:
    exactsol = pickle.load(f)

step = int( exactsol[-1].getNe() / mesh.getNe() )
sol = exactsol[0][::step] # sampling
# ===========================================

from timeStepping import CF_Asch
import Initdata as Initd
bc = Initd.Source('NlinAdv0007',data.Re)

nsteps = 10
timestepping = CF_Asch(data, gll, mesh, bc, nsteps)
timestepping.stepBystepReport = False 
N = 4
err = np.zeros((N,))
dt = np.zeros((N,))

for n in range(N):
    timestepping.evolve()
    err[n] = np.linalg.norm(timestepping.sol - sol, np.inf) 
    dt[n] = timestepping.dt
    print("l2-Error obtain with %i steps = %f." %(nsteps,err[n]))
    nsteps *= 2
    timestepping.nsteps = nsteps
    timestepping.nsteps2do = nsteps
    timestepping.dt = data.dT / nsteps

normsol = np.linalg.norm(sol, np.inf) 
err /= normsol  # relative error
# ===========================================
from matplotlib import pyplot as plt

h1, = plt.plot(mesh.nodes,timestepping.sol)
h2, = plt.plot(exactsol[-1].nodes,exactsol[0])
plt.xlabel('x')
plt.legend([h1,h2],['SLRK3','ref'])
plt.show()

h1, = plt.loglog(dt, err,'.-')
h2, = plt.loglog(dt, dt**3,'--')
plt.xlabel('h')
plt.ylabel('Error (L2)')
plt.legend([h1,h2],['SLRK3','h^3'])
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

print('slopes =', slopes)

with open('temporal1.npy','wb') as f:
    np.save(f, err)
    np.save(f, dt)

