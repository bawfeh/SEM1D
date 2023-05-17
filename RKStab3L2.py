# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
from Methods import DIRK_CF
from scipy.optimize import bisect, fsolve

data = InputFile.Parameters() # default input data
data.cfmethod = data.cfmethod.fromkeys(data.cfmethod,False)
data.cfmethodName = 'CF343'
data.cfmethod[data.cfmethodName ] = True

cf = DIRK_CF( np, data.cfmethod )

gll = GaussMesh.GaussDistribution(40)

def rho(v, i=1):
    return 1 - cf.AI[i,i] * v

def y1(v, w, real=True):
    return 1./rho(v) if real else cf.AE[1,0]*w/rho(v)

def rhoy1(v, w, i, real=True):
    if i < 2: return 0.
    return cf.AE[i,1]*w*y1(v,w,False) - cf.AI[i,1]*v*y1(v,w) \
            if real else \
            -cf.AI[i,1]*v*y1(v,w,False) - cf.AE[i,1]*w*y1(v,w)

def y2(v, w, real=True):
    return (1 - rhoy1(v,w,2)) / rho(v,2) \
            if real else \
            (cf.AE[2,0]*w - rhoy1(v,w,2,False)) / rho(v,2)

def rhoy2(v, w, i, real=True):
    if i < 3: return 0.
    return cf.AE[i,2]*w*y2(v,w,False) - cf.AI[i,2]*v*y2(v,w) \
            if real else \
            -cf.AI[i,2]*v*y2(v,w,False) - cf.AE[i,2]*w*y2(v,w)

def y3(v, w, real=True):
    return (1 - rhoy1(v,w,3) - rhoy2(v,w,3)) / rho(v,3) \
            if real else \
            (cf.AE[3,0]*w - rhoy1(v,w,3,False) - rhoy2(v,w,3,False)) / rho(v,3)

def R(v, w):
    return 1 + v * (cf.b[1]*y1(v,w) + cf.b[2]*y2(v,w) + cf.b[3]*y3(v,w)) \
            - w * (cf.b[1]*y1(v,w,False) + cf.b[2]*y2(v,w,False) + cf.b[3]*y3(v,w,False))

def I(v, w):
    return v * (cf.b[1]*y1(v,w,False) + cf.b[2]*y2(v,w,False) + cf.b[3]*y3(v,w,False)) \
            + w * (cf.b[1]*y1(v,w) + cf.b[2]*y2(v,w) + cf.b[3]*y3(v,w))

def mod(v, w): return I(v,w)**2 + R(v,w)**2


w0 = -5.0; w1 = 5.0
# w = w0 + ((w1-w0)/2.) * (1+gll.gllnodes)
w = np.concatenate((w0 + ((0.-w0)/2.) * (1+gll.gllnodes), \
        0. + ((w1-0.)/2.) * (1+gll.gllnodes)))
w = np.unique(w)
v = np.zeros(np.shape(w))
v0 = -8.; v1 = 4.
val = np.zeros(np.shape(w), dtype=bool)
for i in range(len(w)):
    f = lambda v : R(v,w[i])**2 + I(v,w[i])**2 - 1.
    fab = f(v0) * f(v1)
    if fab > 0: print('>>>', w[i]); continue
    x0, r = bisect(f, v0, v1, full_output=True, xtol=1e-3)
    if not(r.converged): print('///', w[i]); continue
    sol = fsolve(f, x0)#, xtol=1e-14)
    if np.isclose(f(sol[0]),0.):
        val[i] = True; v[i] = sol[0]
    else: print('&&&', w[i])

O =  np.vstack((v[val],w[val])).T
print('Near the origin', O.shape)

# np.savetxt('temp0.dat', O, fmt=['%+14.12f','%+9.5f'])
np.savetxt('temp0.dat', O)
# ==========


