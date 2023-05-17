# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh

filename = '/Users/bawfeh78/Documents/Cplusplus/SLDG/input5.txt'

data = InputFile.Parameters() # default input data

flag = InputFile.InputFile(filename).updateData(data) # with file inputs

# if flag : print('\nData successfully updated!\n')
# else: print('\nWARNING: Incorrect data provided in '+ filename)

from Methods import RKDG

method = RKDG(np, data.orderRK)

def y(t): return np.exp(t)
def f(t,y): return y

# solution to y' = y, y(0) = 1.
N = 5
err = np.zeros((N,))
dt = np.zeros((N,))
nsteps = 10
nstages = len(method.alpha)
T = 1.
Y = [0.]*nstages
F = [0.]*nstages

for n in range(N):
    y0 = 1.
    dt[n] = T/nsteps
    for _ in range(nsteps):
        for i in range(nstages):
            Y[i] = y0
            F[i] = y0
            y0 = method.alpha[i].dot(Y) + dt[n] * method.beta[i].dot(F)
    err[n] = abs(y(T) - y0)
    if n > 0:
        order = (np.log(err[n]) - np.log(err[n-1])) / (np.log(dt[n]) - np.log(dt[n-1]))
        print('Using %d steps, order = %f' %(nsteps, order))
    nsteps *= 2
print("errors =", err)
# ===========================================
from matplotlib import pyplot as plt

h1, = plt.loglog(dt, err,'.-')
h2, = plt.loglog(dt, dt**data.orderRK,'--')
plt.xlabel('h')
plt.ylabel('Error (L1)')
plt.legend([h1,h2],['RKDG','h^%i'%data.orderRK])
plt.show()
