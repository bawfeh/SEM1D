# Numpy tutorial: basics
# import sys, os
import numpy as np

ext = 'eps'
source = '/Users/bawfeh78/Documents/PYTHON/SEM1D/Data/'
dest = ''

# ============================================
# with open(source+'varyViscosity1.npy', 'rb') as f:
#     errs = np.load(f)
#     nu = np.load(f)
# 
# from matplotlib import pyplot as plt
# 
# fig, ax = plt.subplots()
# h1, = ax.loglog(nu, errs[:,0],'.-', lw=2, ms=10)
# h2, = ax.loglog(nu, errs[:,1],'x-', lw=2, ms=6)
# h3, = ax.loglog(nu, errs[:,2],'d-', lw=2, ms=6)
# plt.xlabel(r'$\nu$',fontsize='large',fontweight='semibold')
# plt.ylabel('Rel. Error (L2)',fontsize='large',fontweight='semibold')
# # ax.set(ylim=[3e-3,1.1e-0])
# ax.tick_params(length=6, width=2, labelsize=14, direction='in')
# ax.tick_params(which='minor', length=4, direction='in')
# legend_props = {'weight':'bold','size':12}
# plt.legend([h1,h2,h3],['SLRK3','SLRK3L','SLRK3L(2)'], frameon=False,prop=legend_props,loc='best')
# plt.savefig('errors3.eps', bbox_inches = 'tight')
# plt.show()
# ============================================
# 
# with open(source+'varyViscosity1.npy', 'rb') as f:
#     errs = np.load(f)
#     nu = np.load(f)
# 
# with open(source+'varyViscCF233CK.npy', 'rb') as f:
#     err = np.load(f)
#     nu1 = np.load(f)
# 
# from matplotlib import pyplot as plt
# 
# fig, ax = plt.subplots()
# h1, = ax.loglog(nu, errs[:,0],'.-', lw=2, ms=10)
# h2, = ax.loglog(nu1, err,'x-', lw=2, ms=6)
# plt.xlabel(r'$\nu$',fontsize='large',fontweight='semibold')
# plt.ylabel('Rel. Error (L2)',fontsize='large',fontweight='semibold')
# # ax.set(ylim=[3e-3,1.1e-0])
# ax.tick_params(length=6, width=2, labelsize=14, direction='in')
# ax.tick_params(which='minor', length=4, direction='in')
# legend_props = {'weight':'bold','size':12}
# plt.legend([h1,h2],['SLRK3','SLRK3[CK]'], frameon=False,prop=legend_props,loc='best')
# plt.savefig('errors4.eps', bbox_inches = 'tight')
# plt.show()
# # ============================================

errs = []
with open(source+'varyViscosityIMEX233.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscosityIMEX443.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscosityIMEX343.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL233.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL443.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL343.npy', 'rb') as f:
    err = np.load(f)
    nu = np.load(f)
    errs.append(err)

from matplotlib import pyplot as plt

fig, ax = plt.subplots()
h1, = ax.loglog(nu, errs[0],'b.--', lw=2, ms=10)
h2, = ax.loglog(nu, errs[1],'rx--', lw=2, ms=6)
h3, = ax.loglog(nu, errs[2],'md--', lw=2, ms=6)
h4, = ax.loglog(nu, errs[3],'b.-', lw=2, ms=10)
h5, = ax.loglog(nu, errs[4],'rx-', lw=2, ms=6)
h6, = ax.loglog(nu, errs[5],'md-', lw=2, ms=6)
plt.xlabel(r'$\mathbf{\nu}$',fontsize=18,fontweight='bold')
plt.ylabel(r'Rel.Error (L$_2$)',fontsize=14,fontweight='semibold')
# ax.set(ylim=[3e-3,1.1e-0])
ax.tick_params(length=6, width=2, labelsize=14, direction='in')
ax.tick_params(which='minor', length=4, direction='in')
legend_props = {'weight':'bold','size':14}
plt.legend([h1,h2,h3,h4,h5,h6],['IMEX3','IMEX3L','IMEX3L(2)','SLRK3','SLRK3L','SLRK3L(2)'], frameon=False,prop=legend_props,loc='best')
plt.savefig(dest+'varyViscosityALL.'+ext, bbox_inches = 'tight')
plt.show()
# ============================================

errs = []
with open(source+'varyViscositySL233CK.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL233.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL443CKV.npy', 'rb') as f:
    err = np.load(f)
    errs.append(err)
with open(source+'varyViscositySL443.npy', 'rb') as f:
    err = np.load(f)
    nu = np.load(f)
    errs.append(err)

from matplotlib import pyplot as plt

fig, ax = plt.subplots()
h1, = ax.loglog(nu, errs[0],'b.--', lw=2, ms=10)
h2, = ax.loglog(nu, errs[1],'b.-', lw=2, ms=10)
h3, = ax.loglog(nu, errs[2],'rx--', lw=2, ms=6)
h4, = ax.loglog(nu, errs[3],'rx-', lw=2, ms=6)
plt.xlabel(r'$\mathbf{\nu}$',fontsize=18,fontweight='bold')
plt.ylabel(r'Rel.Error (L$_2$)',fontsize=14,fontweight='semibold')
# ax.set(ylim=[3e-3,1.1e-0])
ax.tick_params(length=6, width=2, labelsize=14, direction='in')
ax.tick_params(which='minor', length=4, direction='in')
legend_props = {'weight':'bold','size':14}
plt.legend([h1,h2,h3,h4],['SLRK3[CK]','SLRK3','SLRK3L[CKV]','SLRK3L'], frameon=False,prop=legend_props,loc='best')
plt.savefig(dest+'varyViscositySLOpt.'+ext, bbox_inches = 'tight')
plt.show()
# ============================================
