# Numpy tutorial: basics
# import sys, os
import numpy as np
from scipy.optimize import least_squares
# import sympy as sym

## --- DIRK(2,3,3) ---
# gma = (3+np.sqrt(3))/6
# A = np.array([[0,0,0],[0,gma,0],[0,1-2*gma,gma]])
# B = np.array([[0,0,0],[gma,0,0],[gma-1,2*(1-gma),0]])
# c = np.array([0,gma,1-gma])

## --- DIRK(3,3,3) ---
# bta = np.sqrt(3)/3
# A = np.array([[0,0,0], [-bta/2,(1+bta)/2,0], [(3+5*bta)/2,-(1+3*bta),(1+bta)/2]])
# B = np.array([[0,0,0], [.5,0,0], [-1,2,0]])
# c = np.array([0,.5,1])

## --- DIRK(4,4,3) ---
# A = np.array([[0,0,0,0,0],[0,.5,0,0,0],[0,1./6,.5,0,0],[0,-.5,.5,.5,0],[0,1.5,-1.5,.5,.5]])
# B = np.array([[0,0,0,0,0],[.5,0,0,0,0],[11./18,1./18,0,0,0],[5./6,-5./6,.5,0,0],[1./4,7./4,3./4,-7./4,0]])
# c = np.array([0,.5,2./3,.5,1])

## --- DIRK(3,4,3) ---
gma = 0.4358665215
dta = 0.5529291479
b1 = -3*gma**2/2+4*gma-1./4
a31 = (15./4-15*gma+21*gma**2/4)*dta - 7./2+13*gma-9*gma**2/2
A = np.array([[0.,0.,0.,0.],[0.,gma,0.,0.],[0.,(1-gma)/2,gma,0.],[0.,b1,1-gma-b1,gma]])
B = np.array([[0.,0.,0.,0.],[gma,0.,0.,0.],[a31,(1+gma)/2-a31,0.,0.],[1-2*dta,dta,dta,0.]]) 
c = np.array([0.,gma,(1+gma)/2,1.])

# print('A =', A)
# print('B =', B)
# print('c =', c)

# def F(x):
#     y = np.zeros((5,))
#     ones = np.ones(np.shape(c))
#     x[-1] = (1./3 - x[:-1].dot(ones[:-1]-2*c[:-1]))/(1-2*c[-1])
#     y[0] = -7 * x.dot(c) * np.sum(x) + 2 * np.sum(x)**2 + 7 * x.dot(c) - 2 * np.sum(x) - 1./2
#     y[1] =  2 * x.dot((ones-c)*c) * np.sum(x) - x.dot((ones+c)*c) + 1./2 * np.sum(x) - 1./6
#     y[2] = -2 * x.dot(B.dot(c)) + 1./3 * np.sum(x) - 1./6
#     y[3] = -2 * x.dot(A.dot(c)) + 1./3 * np.sum(x) - 1./6
#     y[4] = -2 * x.dot((c*A).dot(ones)-A.dot(c)) + 1./3 * np.sum(x) - 1./6
#     return y
# 
# x0 = np.zeros(np.shape(c))
# x0[-1] = (1./3 - x0[:-1].dot(np.ones(np.shape(c[:-1]))-2*c[:-1]))/(1-2*c[-1])
# 
# res = least_squares(F,x0,method='lm',max_nfev=8000)
# 
# print(res)

ones = np.ones(np.shape(c))
x = np.zeros(np.shape(c))
def F(xs):
    x[:-1] = xs
    x[-1] = (1./3 - xs.dot(ones[:-1]-2*c[:-1]))/(1-2*c[-1])
    y = np.zeros((5,))
    y[0] = -7 * x.dot(c) * np.sum(x) + 2 * np.sum(x)**2 + 7 * x.dot(c) - 2 * np.sum(x) - 1./2
    y[1] =  2 * x.dot((ones-c)*c) * np.sum(x) - x.dot((ones+c)*c) + 1./2 * np.sum(x) - 1./6
    y[2] = -2 * x.dot(B.dot(c)) + 1./3 * np.sum(x) - 1./6
    y[3] = -2 * x.dot(A.dot(c)) + 1./3 * np.sum(x) - 1./6
    y[4] = -2 * x.dot((c*A).dot(ones)-A.dot(c)) + 1./3 * np.sum(x) - 1./6
    return y

x0 = np.zeros(np.shape(c[:-1]))
# x0 = -np.ones(np.shape(c[:-1]))

res = least_squares(F,x0,max_nfev=8000)

print(res)

x3 = (1./3 - res.x.dot(ones[:-1]-2*c[:-1]))/(1-2*c[-1])

print('xs =', x3)

# ss = res.x.dot(np.ones((3,))-2*c)
# 
# print('constraint = ',ss)

# print('f = ', F(res.x))

# g = sym.Rational(1,2)+sym.sqrt(3)/6
# x1,x2 = sym.symbols('x1, x2')
# c = [0,g,1-g]
# A = sym.Matrix([[0,0,0],[0,g,0],[0,1-2*g,g]])
# B = sym.Matrix([[0,0,0],[g,0,0],[g-1,2*(1-g),0]])
# x = [x1,x2,sym.Rational(1,3)]
# for i in range(2):
#     x[-1] += x[i]*(2*c[i]-1)
# x[-1] /= 1-2*c[-1]
# 
# f = [0]*5
# 
# for k in range(3):
#     f[0] += (7*c[k]-2)*x[k]
#     f[1] += (sym.Rational(1,2)-c[k]*(c[k]+1))*x[k]
#     f[2] += sym.Rational(1,3)*x[k]
#     f[3] += sym.Rational(1,3)*x[k]
#     f[4] += sym.Rational(1,3)*x[k]
#     for j in range(3):
#         f[0] += -7*(c[j]-1)*x[j]*x[k]
#         f[1] += 2*c[k]*(1-c[k])*x[j]*x[k]
#         f[2] += -2*B[j,k]*c[k]*x[j]
#         f[3] += -2*A[j,k]*c[k]*x[j]
#         f[4] += -2*A[j,k]*(c[j]-c[k])*x[j]
# f[0] -= sym.Rational(1,2)
# f[1] -= sym.Rational(1,6)
# f[2] -= sym.Rational(1,6)
# f[3] -= sym.Rational(1,6)
# f[4] -= sym.Rational(1,6)
# 
# J = [[0,0]]*5
# ap = [0]*2
# for m in range(2): ap[m] = -(2*c[m]-1)/(2*c[-1]-1)
# 
# for m in range(2):
#     J[0][m] += 7*(c[m]+ap[m]*c[-1]) - 2*(1+ap[m])
#     J[1][m] += -(1+c[m])*c[m] - ap[m]*(1+c[-1])*c[-1] + sym.Rational(1,2)*(1+ap[m])
#     J[2][m] += sym.Rational(1,3)*(1+ap[m])
#     J[3][m] += sym.Rational(1,3)*(1+ap[m])
#     J[4][m] += sym.Rational(1,3)*(1+ap[m])
#     for k in range(3):
#         J[0][m] += -7*(c[m]+ap[m]*c[-1])*x[k] - 7*(1+ap[m])*c[k]*x[k] + 4*(1+ap[m])*x[k]
#         J[1][m] += 2*(1+ap[m])*c[k]*(1-c[k])*x[k] + c[m]*(1-c[m])*x[k] + ap[m]*c[-1]*(1-c[-1])
#         J[2][m] += -2*(B[m,k]+ap[m]*B[-1,k])*c[k]
#         J[3][m] += -2*(A[m,k]+ap[m]*B[-1,k])*c[k]
#         J[4][m] += -2*A[m,k]*(c[m]-c[k]) - 2*ap[m]*A[-1,k]*(c[-1]-c[k])
# 
# Eq = [0]*2
# 
# for m in range(2):
#     for n in range(5):
#         Eq[m] += J[n][m] * f[n]
# for m in range(2):
#     Eq[m] = sym.simplify(Eq[m])
# 
# solution = sym.solve((Eq[0],Eq[1]),(x1,x2))
# 
# print(solution)
