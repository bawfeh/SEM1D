# Numpy tutorial: basics
# import sys, os
import numpy as np
import InputFile, GaussMesh
from Methods import DIRK_CF
from scipy.optimize import bisect, fsolve

data = InputFile.Parameters() # default input data
data.cfmethod = data.cfmethod.fromkeys(data.cfmethod,False)
data.cfmethodName = 'CF443'
data.cfmethod[data.cfmethodName ] = True

cf = DIRK_CF( np, data.cfmethod )

gll = GaussMesh.GaussDistribution(39)

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

def rhoy3(v, w, i, real=True):
    if i < 4: return 0.
    return cf.AE[i,3]*w*y3(v,w,False) - cf.AI[i,3]*v*y3(v,w) \
            if real else \
            -cf.AI[i,3]*v*y3(v,w,False) - cf.AE[i,3]*w*y3(v,w)

def y4(v, w, real=True):
    return (1 - rhoy1(v,w,4) - rhoy2(v,w,4) - rhoy3(v,w,4)) / rho(v,4) \
            if real else \
            (cf.AE[4,0]*w - rhoy1(v,w,4,False) - rhoy2(v,w,4,False) - rhoy3(v,w,4,False)) / rho(v,4)

b = np.sum(cf.beta,0)
def R(v, w):
    return 1 + v * (cf.b[0] + cf.b[1]*y1(v,w) + cf.b[2]*y2(v,w) + cf.b[3]*y3(v,w) + cf.b[4]*y4(v,w)) \
            - w * (b[1]*y1(v,w,False) + b[2]*y2(v,w,False) + b[3]*y3(v,w,False) + b[4]*y4(v,w,False))

def I(v, w):
    return v * (cf.b[1]*y1(v,w,False) + cf.b[2]*y2(v,w,False) + cf.b[3]*y3(v,w,False) + cf.b[4]*y4(v,w,False)) \
            + w * (b[0] + b[1]*y1(v,w) + b[2]*y2(v,w) + b[3]*y3(v,w) + b[4]*y4(v,w))


w0 = -6.0; w1 = 6.0
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

# u = np.linspace(-8.,-5.0,20)
# x = []; y = []
# for i in range(len(u)):
#     def f(w): return R(u[i],w)**2 + I(u[i],w)**2 - 1.
#     x0 = bisect(f, 0., 8., xtol=1e-6)
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
# u = np.linspace(-5.0,-8.,20)
# x = []; y = []
# for i in range(len(u)):
#     def f(w): return R(u[i],w)**2 + I(u[i],w)**2 - 1.
#     x0 = bisect(f, -8., 0., xtol=1e-6)
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
# w = np.concatenate((np.linspace(-4.5,0,40),np.linspace(0,4.5,40)))
# w = np.unique(w)
# x = []; y = []; x0 = 0.
# for i in range(len(w)):
#     f = lambda u : R(u,w[i])**2 + I(u,w[i])**2 - 1.
#     x0, r = bisect(f, -8., 4., full_output=True, xtol=1e-6)
#     if r.converged:
#         sol = fsolve(f, x0)#, xtol=1e-14)
#         x.append(sol[0]); y.append(w[i])
#     else: print('>>>', w[i])
# 
# O = np.flipud(np.array([x, y]).T)
# print('Near the origin', O.shape)
# 
# # print('%8s %8s' % ('x','y'))
# for d in O:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# print('='*30)
# # ==========

# # w = np.concatenate((np.linspace(-6.5,0,40),np.linspace(0,6.5,40)))
# # w = np.unique(w)
# w0 = -6.0; w1 = 6.0
# w = w0 + ((w1-w0)/2.) * (1+gll.gllnodes)
# x = []; y = []; x0 = 0.
# for i in range(len(w)):
#     f = lambda u : R(u,w[i])**2 + I(u,w[i])**2 - 1.
#     fa = f(-8.); fb = f(4.)
#     if fa*fb > 0: continue
#     x0, r = bisect(f, -8., 4., full_output=True, xtol=1e-6)
#     if r.converged:
#         sol = fsolve(f, x0)#, xtol=1e-14)
#         x.append(sol[0]); y.append(w[i])
#     else: print('>>>', w[i])
# 
# O = np.flipud(np.array([x, y]).T)
# print('Near the origin', O.shape)
# 
# # print('%8s %8s' % ('x','y'))
# for d in O:
#     print('%8.5f %8.5f' % (d[0], d[1]))
# print('='*30)
# ==========
# # 
# # u = np.concatenate((np.linspace(-8.0,0,50),np.linspace(0,1,50)))
# # u = np.concatenate((u,np.linspace(1,5,50)))
# # # u = np.concatenate((np.linspace(0.92500,1,50),np.linspace(1,2,5)))
# # u = np.unique(u)
# # x = []; y = []
# # for s in u:
# #     def f(v): return R(v,0.)**2 + I(v,0.)**2 - 1.
# #     if f(s) <= 0.:
# #         x.append(s); y.append(0.)
# #     else: print('>>>>',s)
# # 
# # print('='*30)
# # print(R(0.95833,0.)**2 + I(0.95833,0.)**2 - 1.)
# # A = np.array([x, y]).T
# # print('Implicit RK', A.shape)
# # 
# # # print('%8s %8s' % ('x','y'))
# # for d in A:
# #     print('%8.5f %8.5f' % (d[0], d[1]))
# # 
# # # v0 = 0.1; w0 = 0.
# # # for i in range(10):
# # #     w0 += 0.2
# # #     print( R(v0,w0)**2 + I(v0,w0)**2 - 1 )
# # 
# # # print(R(U[:,0],U[:,1])**2 + I(U[:,0],U[:,1])**2 - 1)
# # 
# # # with open('temp1.npy', 'wb') as mf:
# # #     np.save(mf, U)
# # # 
