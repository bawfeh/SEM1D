# Numpy tutorial: basics
# import sys, os
import numpy as np

ext = 'eps'
source = '/Users/bawfeh78/Documents/PYTHON/SEM1D/Data/'
dest = ''

# # ============================================
# with open(source+'temporalAccuracy.npy', 'rb') as f:
#     errs = np.load(f)
#     dt = np.load(f)
# 
# from matplotlib import pyplot as plt
# 
# fig, ax = plt.subplots()
# h1, = ax.loglog(dt, errs[:,0],'.-', lw=2, ms=10)
# h2, = ax.loglog(dt, errs[:,1],'x-', lw=2, ms=6)
# h3, = ax.loglog(dt, errs[:,2],'d-', lw=2, ms=6)
# h4, = ax.loglog(dt, dt**3,'purple')
# plt.xlabel('h',fontsize='large',fontweight='semibold')
# plt.ylabel('Rel. Error (L2)',fontsize='large',fontweight='semibold')
# # ax.set(ylim=[1e-7,1.1e-1])
# ax.set(xlim=[1e-3,2.5e-1])
# ax.tick_params(length=6, width=2, labelsize=14, direction='in')
# ax.tick_params(which='minor', length=4, direction='in')
# legend_props = {'weight':'bold','size':12}
# plt.legend([h1,h2,h3,h4],['SLRK3','SLRK3L','SLRK3L(2)','h^3'], frameon=False,prop=legend_props,loc='best')
# plt.savefig('errors2.eps', bbox_inches = 'tight')
# plt.show()
# # ============================================

errs = []
with open(source+'temporalAccuracySL233.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'temporalAccuracySL443.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'temporalAccuracySL343.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
    dt = np.load(f)

from matplotlib import pyplot as plt

fig, ax = plt.subplots()
h1, = ax.loglog(dt, errs[0],'b.-', lw=2, ms=10)
h2, = ax.loglog(dt, errs[1],'rx-', lw=2, ms=6)
h3, = ax.loglog(dt, errs[2],'md-', lw=2, ms=6)
h4, = ax.loglog(dt, dt**3,'k-')
plt.xlabel(r'$\mathbf{h}$',fontsize=18,fontweight='bold')
plt.ylabel(r'Rel.Error (L$_2$)',fontsize=14,fontweight='semibold')
# ax.set(ylim=[1e-7,1.1e-1])
# ax.set(xlim=[1e-3,2.5e-1])
ax.tick_params(length=6, width=2, labelsize=14, direction='in')
ax.tick_params(which='minor', length=4, direction='in')
legend_props = {'weight':'bold','size':12}
plt.legend([h1,h2,h3,h4],['SLRK3','SLRK3L','SLRK3L(2)',r'h$^3$'], frameon=False,prop=legend_props,loc='lower right')
plt.savefig(dest+'temporalAccuracySL.'+ext, bbox_inches = 'tight')
plt.show()
# ============================================

errs = []
with open(source+'temporalAccuracySL233CK.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'temporalAccuracySL233.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'temporalAccuracySL443CKV.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'temporalAccuracySL443.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
    dt = np.load(f)

from matplotlib import pyplot as plt

fig, ax = plt.subplots()
h1, = ax.loglog(dt, errs[0],'b.--', lw=2, ms=10)
h2, = ax.loglog(dt, errs[1],'b.-', lw=2, ms=10)
h3, = ax.loglog(dt, errs[2],'rx--', lw=2, ms=6)
h4, = ax.loglog(dt, errs[3],'rx-', lw=2, ms=6)
h5, = ax.loglog(dt, dt**3,'k-')
plt.xlabel(r'$\mathbf{h}$',fontsize=18,fontweight='bold')
plt.ylabel(r'Rel.Error (L$_2$)',fontsize=14,fontweight='semibold')
# ax.set(ylim=[1e-7,1.1e-1])
# ax.set(xlim=[1e-3,2.5e-1])
ax.tick_params(length=6, width=2, labelsize=14, direction='in')
ax.tick_params(which='minor', length=4, direction='in')
legend_props = {'weight':'bold','size':12}
plt.legend([h1,h2,h3,h4,h5],['SLRK3[CK]','SLRK3','SLRK3L[CKV]','SLRK3L',r'h$^3$'], frameon=False,prop=legend_props,loc='lower right')
plt.savefig(dest+'temporalAccuracySLOpt.'+ext, bbox_inches = 'tight')
plt.show()
# ============================================

# slopes = np.log(errs[1:,:])-np.log(errs[:-1,:])
# hh = np.log(dt[1:])-np.log(dt[:-1])
# for k in range(3):
#     slopes[:,k] /= hh
# 
# print(slopes)
