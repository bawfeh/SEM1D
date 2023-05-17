# Numpy tutorial: basics
# import sys, os
import numpy as np

# folder = '/Users/bawfeh78/Documents/Research/Manuscripts2/Figures/'
folder = '/Users/bawfeh78/Documents/Research/Manuscripts3/Figures/'
outputf = './LinearTimeRKDG/'

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, NullFormatter )
from matplotlib import ticker

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

with open(outputf+'LinearTimeRKDG1pConvCr2p5.npy', 'rb') as f:
    errs = np.load(f)
    qq = np.load(f)

for ext in ['.pdf','.eps']:
    fig, ax = plt.subplots()
    h1, = ax.semilogy(qq, errs[0], 'b.-', lw=2, ms=8)
    h2, = ax.semilogy(qq, errs[1], 'rx-', lw=2, ms=6)
    h3, = ax.semilogy(qq, errs[2], 'ms-', lw=2, ms=6)
    plt.xlabel('q',fontsize=16,fontweight='semibold')
    plt.ylabel(r'Error (L$_2$)', fontsize=18,fontweight='semibold')
    # ax.set(xlim=[0,2.1])
    ax.set(ylim=[1e-11,1e-1])
    locmin = ticker.LogLocator(base=10, numticks=14)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.tick_params(length=6, width=2, labelsize=14, direction='in')
    ax.tick_params(which='minor', length=4, direction='in')
    legend_props = {'weight':'bold','size':14}
    plt.legend([h1,h2,h3],['M = 20','M = 40', 'M = 80'], frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsLinearTimeRKDG1hpConv'+ext, bbox_inches = 'tight')
    plt.show()

# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
with open(outputf+'LinearTimeRKDG1hConvCr10.npy', 'rb') as f:
    errs = np.load(f)
    M = np.load(f)
    dt = np.load(f)

for ext in ['.pdf','.eps']:
    fig, ax = plt.subplots()
    h1, = ax.loglog(2./M, errs[0], 'b.-', lw=2, ms=8)
    h2, = ax.loglog(2./M, errs[1], 'rx-', lw=2, ms=6)
    h3, = ax.loglog(2./M, errs[2], 'ms-', lw=2, ms=6)
    plt.xlabel('h',fontsize=16,fontweight='semibold')
    plt.ylabel(r'Error (L$_2$)', fontsize=18,fontweight='semibold')
    ax.set(ylim=[1e-9,1e-1])
    ax.tick_params(length=6, width=2, labelsize=14, direction='in')
    ax.tick_params(which='minor', length=4, direction='in')
    legend_props = {'weight':'bold','size':14}
    plt.legend([h1,h2,h3],[r'$\mathbf{q = 1}$',r'$\mathbf{q = 2}$',r'$\mathbf{q = 3}$'],\
            frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsLinearTimeRKDG1hpCr10Conv'+ext, bbox_inches = 'tight')
    plt.show()

# ===========================================

errs *= dt
slopes = (np.log(errs[:,1:])-np.log(errs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
print('slopes =',slopes)

