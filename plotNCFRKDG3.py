# Numpy tutorial: basics
# import sys, os
import numpy as np
# ext = '.eps'
ext = '.pdf'
# folder = '/Users/bawfeh78/Documents/Research/Manuscripts2/Figures/'
folder = '/Users/bawfeh78/Documents/Research/Manuscripts3/Figures/'
outputf = './NCFRKDG/'

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

# ============================================

with open(outputf+'NCFRKDG3.npy', 'rb') as f:
    solA = np.load(f)
    x = np.load(f)
    t = np.load(f)
    exactsol = np.load(f)

with open(outputf+'NCFRKDG3SL.npy', 'rb') as f:
    sol = np.load(f)

for ext in ['.pdf','.eps']:
    figname = ['sol1NCFRKDG3xAB'+ext,\
            'sol2NCFRKDG3xAB'+ext,\
            'sol3NCFRKDG3xAB'+ext,\
            'sol4NCFRKDG3xAB'+ext]
    
    for k in range(4):
        print('t =', t[k+1])
        fig, ax = plt.subplots()
        h1, = ax.plot(x, exactsol[k+1], lw=2, ms=6)
        h2, = ax.plot(x, sol[k+1],'m--', lw=2, ms=6)
        h3, = ax.plot(x, solA[k+1],'r', lw=1, ms=6)
        ax.set_xlabel('x',fontsize=16,fontweight='semibold')
        ax.set_ylabel('u',fontsize=16,fontweight='semibold')
        ax.yaxis.set_major_locator(MultipleLocator(.2))
        ax.yaxis.set_minor_locator(MultipleLocator(.1))
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.xaxis.set_minor_locator(MultipleLocator(.1))
        ax.tick_params(length=6, width=2, labelsize=14, direction='in')
        ax.tick_params(which='minor', length=4, direction='in')
        legend_props = {'weight':'bold','size':14}
        ax.legend([h1,h2,h3],['Ref.sol.','SLRKDG3', 'RKDG3'], frameon=False,prop=legend_props,loc='best')
#         plt.savefig(folder+figname[k], bbox_inches = 'tight')
    plt.show()
# ============================================

# with open(outputf+'NCFRKDG3SL.npy', 'rb') as f:
#     sol = np.load(f)
# 
# with open(outputf+'NCFRKDG3SL2.npy', 'rb') as f:
#     sol1 = np.load(f)
# 
# with open(outputf+'NCFRKDG3SL1.npy', 'rb') as f:
#     sol2 = np.load(f)
# 
# for ext in ['.pdf','.eps']:
#     figname = ['sol1NCFRKDG3xB'+ext,\
#             'sol2NCFRKDG3xB'+ext,\
#             'sol3NCFRKDG3xB'+ext,\
#             'sol4NCFRKDG3xB'+ext]
#     
#     for k in range(4):
#         print('t =', t[k+1])
#         fig, ax = plt.subplots()
#         h1, = ax.plot(x, exactsol[k+1], lw=2, ms=6)
#         h2, = ax.plot(x, sol[k+1],'m--', lw=2, ms=6)
#         h3, = ax.plot(x, sol1[k+1],'r', lw=1, ms=6)
#         h4, = ax.plot(x, sol2[k+1],'k-.', lw=1, ms=6)
#         plt.xlabel('x',fontsize=16,fontweight='semibold')
#         plt.ylabel('u',fontsize=16,fontweight='semibold')
#         ax.yaxis.set_major_locator(MultipleLocator(.2))
#         ax.yaxis.set_minor_locator(MultipleLocator(.1))
#         ax.xaxis.set_major_locator(MultipleLocator(.5))
#         ax.xaxis.set_minor_locator(MultipleLocator(.1))
#         ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#         ax.tick_params(which='minor', length=4, direction='in')
#         legend_props = {'weight':'bold','size':14}
#         plt.legend([h1,h2,h3,h4],['Ref.sol.','case 1', 'case 2', 'case 3'], frameon=False,prop=legend_props,loc='best')
# #         plt.savefig(folder+figname[k], bbox_inches = 'tight')
#     plt.show()
# # ============================================
