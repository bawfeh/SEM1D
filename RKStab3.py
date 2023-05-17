# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
from Methods import DIRK_CF
from scipy.optimize import bisect, fsolve

data = InputFile.Parameters() # default input data
data.cfmethod = data.cfmethod.fromkeys(data.cfmethod,False)
data.cfmethodName = 'CF233'
data.cfmethod[data.cfmethodName ] = True

cf = DIRK_CF( np, data.cfmethod )

gll = GaussMesh.GaussDistribution(24)

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

def R(v, w):
    return 1 + v * (cf.b[1]*y1(v,w) + cf.b[2]*y2(v,w)) \
            - w * (cf.b[1]*y1(v,w,False) + cf.b[2]*y2(v,w,False))

def I(v, w):
    return v * (cf.b[1]*y1(v,w,False) + cf.b[2]*y2(v,w,False)) \
            + w * (cf.b[1]*y1(v,w) + cf.b[2]*y2(v,w))


v = np.linspace(-8.,-2.0,50)
w0 = 0.; w1 = 8.; x0 = 0.
w = np.zeros(np.shape(v))
val = np.zeros(np.shape(v), dtype=bool)
for i in range(len(v)):
    def f(w): return R(v[i],w)**2 + I(v[i],w)**2 - 1.
    fab = f(w0) * f(w1)
    if fab > 0: print('>>>', w[i]); continue
    x0, r = bisect(f, w0, w1, full_output=True, xtol=1e-3)
    if not(r.converged): print('>>>', v[i]); continue
    sol = fsolve(f, x0)#, xtol=1e-14)
    if np.isclose(f(sol[0]),0.):
        val[i] = True; w[i] = sol[0]
    else: print('>>>', v[i])

U =  np.vstack((v[val],w[val])).T
print('Upper bound', U.shape)
print('='*30)
# ==========

w0 = -2.5; w1 = 2.5
w = np.concatenate((w0 + ((0.-w0)/2.) * (1+gll.gllnodes), \
        0. + ((w1-0.)/2.) * (1+gll.gllnodes)))
w = np.unique(w)
v = np.zeros(np.shape(w))
v0 = -2.; v1 = 1.; x0 = 0.
val = np.zeros(np.shape(w), dtype=bool)
for i in range(len(w)):
    f = lambda v : R(v,w[i])**2 + I(v,w[i])**2 - 1.
    sol = fsolve(f, x0)#, xtol=1e-14)
    if np.isclose(f(sol[0]),0.):
        val[i] = True; v[i] = sol[0]
    else: print('>>>', w[i])

O =  np.flipud(np.vstack((v[val],w[val])).T)
print('Near the origin', O.shape)
print('='*30)
# ==========

v = np.linspace(-2.0,-8.,50)
w0 = -8.; w1 = 0.; x0 = 0.
w = np.zeros(np.shape(v))
val = np.zeros(np.shape(v), dtype=bool)
for i in range(len(v)):
    def f(w): return R(v[i],w)**2 + I(v[i],w)**2 - 1.
    fab = f(w0) * f(w1)
    if fab > 0: print('>>>', w[i]); continue
    x0, r = bisect(f, w0, w1, full_output=True, xtol=1e-3)
    if not(r.converged): print('>>>', v[i]); continue
    sol = fsolve(f, x0)#, xtol=1e-14)
    if np.isclose(f(sol[0]),0.):
        val[i] = True; w[i] = sol[0]
    else: print('>>>', v[i])

D =  np.vstack((v[val],w[val])).T
print('Lower bound', D.shape)
print('='*30)
# ==========


np.savetxt('temp0.dat', U)
with open('temp0.dat','ab') as mfile:
    np.savetxt(mfile, O)
    np.savetxt(mfile, D)


# 
# u = np.linspace(-8.,-2.0,50)
# x = []; y = []
# for i in range(len(u)):
#     def f(w): return R(u[i],w)**2 + I(u[i],w)**2 - 1.
#     x0, r = bisect(f, 0., 8., full_output=True)
#     sol = fsolve(f, x0)#, xtol=1e-14)
#     if np.isclose(f(sol[0]),0.):
#         x.append(u[i]); y.append(sol[0])
#     else: print('>>>', u[i])
# 
# U = np.array([x, y]).T
# print('Upper bound', U.shape)
# 
# # print('%8s %8s' % ('x','y'))
# for d in U:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# print('='*30)
# # ==========
# 
# u = np.linspace(-2.0,-8.,50)
# x = []; y = []
# for i in range(len(u)):
#     def f(w): return R(u[i],w)**2 + I(u[i],w)**2 - 1.
#     x0, r = bisect(f, -8., 0., full_output=True, xtol=1e-6)
#     sol = fsolve(f, x0)#, xtol=1e-14)
#     if np.isclose(f(sol[0]),0.):
#         x.append(u[i]); y.append(sol[0])
#     else: print('>>>', u[i])
# 
# D = np.array([x, y]).T
# print('Lower bound', D.shape)
# 
# # print('%8s %8s' % ('x','y'))
# for d in D:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# print('='*30)
# # ==========
# 
# w = np.concatenate((np.linspace(-2.5,0,25),np.linspace(0,2.5,25)))
# w = np.unique(w)
# x = []; y = []
# for i in range(len(w)):
#     def f(u): return R(u,w[i])**2 + I(u,w[i])**2 - 1.
#     sol = fsolve(f, 0.,xtol=1e-8)
#     if np.isclose(f(sol[0]),0.):
#         x.append(sol[0]); y.append(w[i])
#     else: print('>>>', w[i])
# # 
# O = np.flipud(np.array([x, y]).T)
# print('Near the origin', O.shape)
# 
# # print('%8s %8s' % ('x','y'))
# for d in O:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# print('='*30)
# # ==========
# print('All U < O < L')
# for d in U:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# for d in O:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# for d in D:
#     print('%8.5f %8.5f' % (d[0], d[1]))
