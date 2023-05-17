# Numpy tutorial: basics
# import sys, os
import numpy as np
# ext = '.eps'
ext = '.pdf'
folder = '/Users/bawfeh78/Documents/Research/Manuscripts2/Figures/'
# folder = '/Users/bawfeh78/Documents/Research/Manuscripts3/Figures/'
outputf = './BurgersGRKDG/'

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

# # ============================================
# 
# with open(outputf+'BurgersGRKDG3.npy', 'rb') as f:
#     sol = np.load(f)
#     x = np.load(f)
#     t = np.load(f)
#     exactsol = np.load(f)
#     errsA = np.load(f)
# 
# with open(outputf+'BurgersGRKDG3SL.npy', 'rb') as f:
#     solB = np.load(f)
#     x = np.load(f)
#     t = np.load(f)
#     exactsol = np.load(f)
#     errsB = np.load(f)
# 
# for ext in ['.pdf','.eps']:
#     figname = ['sol1BurgersGRKDG3xAB'+ext,\
#             'sol2BurgersGRKDG3xAB'+ext,\
#             'sol3BurgersGRKDG3xAB'+ext,\
#             'sol4BurgersGRKDG3xAB'+ext]
#     loc = ['lower right', 'lower right', 'best', 'best']
#     
#     for k in range(4):
#         print('t =', t[5*(k+1)])
#         fig, ax = plt.subplots()
#         h1, = ax.plot(x, exactsol[5*(k+1)], lw=2, ms=6)
#         h2, = ax.plot(x, sol[5*(k+1)],'r', lw=1, ms=6)
#         h3, = ax.plot(x, solB[5*(k+1)],'m--', lw=2, ms=6)
#         plt.xlabel('x',fontsize=16,fontweight='semibold')
#         plt.ylabel('u',fontsize=16,fontweight='semibold')
#         ax.set(ylim=[.2,.8])
#         ax.yaxis.set_major_locator(MultipleLocator(.1))
#         ax.yaxis.set_minor_locator(MultipleLocator(.05))
#         ax.xaxis.set_major_locator(MultipleLocator(.5))
#         ax.xaxis.set_minor_locator(MultipleLocator(.1))
#         ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#         ax.tick_params(which='minor', length=4, direction='in')
#         legend_props = {'weight':'bold','size':14}
#         plt.legend([h1,h2,h3],['Ref.sol.','RKKDG3','SLRKKDG3'], frameon=False,prop=legend_props,loc=loc[k])
# #         plt.savefig(folder+figname[k], bbox_inches = 'tight')
#     plt.show()
# # ============================================
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.semilogy(t[1:], errsA[1:], 'r.-', lw=2, ms=8)
#     h2, = ax.semilogy(t[1:], errsB[1:], 'md-', lw=2, ms=6)
#     plt.xlabel('t',fontsize=16,fontweight='semibold')
#     plt.ylabel('Rel.error (L2)', fontsize=18,fontweight='semibold')
#     ax.set(xlim=[0,2.1])
#     ax.set(ylim=[1e-4,4e-2])
#     ax.xaxis.set_major_locator(MultipleLocator(.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(.1))
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':14}
#     plt.legend([h1,h2],['RKDG','SLRKDG'], frameon=False,prop=legend_props,loc='lower right')
# #     plt.savefig(folder+'errorsBurgersGRKDG3xAB'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# with open(outputf+'BurgersGRKDG3v2.npy', 'rb') as f:
#     sol = np.load(f)
#     x = np.load(f)
#     t = np.load(f)
#     exactsol = np.load(f)
#     errsA = np.load(f)
# 
# with open(outputf+'BurgersGRKDG3SLv2.npy', 'rb') as f:
#     sol = np.load(f)
#     x = np.load(f)
#     t = np.load(f)
#     exactsol = np.load(f)
#     errsB = np.load(f)
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.semilogy(t[1:], errsA[1:], 'r.-', lw=2, ms=8)
#     h2, = ax.semilogy(t[1:], errsB[1:], 'md-', lw=2, ms=6)
#     plt.xlabel('t',fontsize=16,fontweight='semibold')
#     plt.ylabel('Rel.error (L2)', fontsize=18,fontweight='semibold')
#     ax.set(xlim=[0,2.1])
#     ax.set(ylim=[1e-4,4e-2])
#     ax.xaxis.set_major_locator(MultipleLocator(.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(.1))
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':14}
#     plt.legend([h1,h2],['RKDG','SLRKDG'], frameon=False,prop=legend_props,loc='lower right')
# #     plt.savefig(folder+'errorsBurgersGRKDG3xAB2'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# with open(outputf+'BurgersGRKDGRKDG1hConvSL.npy', 'rb') as f:
#     errs = np.load(f)
#     M = np.load(f)
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.loglog(2./M, errs[0], '.-b', lw=2, ms=8)
#     h2, = ax.loglog(2./M, errs[1], 'v-g', lw=2, ms=6)
#     h3, = ax.loglog(2./M, errs[2], 'x-r', lw=2, ms=6)
#     h4, = ax.loglog(2./M, errs[3], 'd-',color='darkorange', lw=2, ms=6)
#     h5, = ax.loglog(2./M, errs[4], 's-m', lw=2, ms=6)
#     plt.xlabel('h',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Rel.error (L$_2$)', fontsize=18,fontweight='semibold')
#     ax.set(ylim=[1e-5,1e-1])
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':12}
#     plt.legend([h1,h2,h3,h4,h5],[r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$'], frameon=False,prop=legend_props,loc='lower right')
#     plt.savefig(folder+'errorsBurgersGRKDG1hpCrConv'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# # ============================================
# with open(outputf+'BurgersGRKDGRKDG3hConvSL.npy', 'rb') as f:
#     errs = np.load(f)
# with open(outputf+'BurgersGRKDGRKDG3hConv.npy', 'rb') as f:
#     errsA = np.load(f)
# 
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.loglog(2./M, errs[0], '.-b', ms=8)
#     h2, = ax.loglog(2./M, errs[1], 'v-g', ms=6)
#     h3, = ax.loglog(2./M, errs[2], 'x-r', ms=6)
#     h4, = ax.loglog(2./M, errs[3], 'd-', color='darkorange', ms=6)
#     h5, = ax.loglog(2./M, errs[4], 's-m', lw=2, ms=6)
#     h6, = ax.loglog(2./M, errsA[0], '.--b', ms=8)
#     h7, = ax.loglog(2./M, errsA[1], 'v--g', ms=6)
#     h8, = ax.loglog(2./M, errsA[2], 'x--r', ms=6)
#     h9, = ax.loglog(2./M, errsA[3], 'd--', color='darkorange', lw=2, ms=6)
#     h10, = ax.loglog(2./M, errsA[4], 's--m', lw=2, ms=6)
#     plt.xlabel('h',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Rel.error (L$_2$)', fontsize=18,fontweight='semibold')
#     ax.set(ylim=[1e-5,1e-1])
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':12}
#     plt.legend([h6,h7,h8,h9,h10,h1,h2,h3,h4,h5],[r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',], frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsBurgersGRKDG3hpCrConvAB'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# # ============================================
# with open(outputf+'BurgersGRKDGRKDG2hConvSL.npy', 'rb') as f:
#     errs = np.load(f)
# with open(outputf+'BurgersGRKDGRKDG2hConv.npy', 'rb') as f:
#     errsA = np.load(f)
# 
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.loglog(2./M, errs[0], '.-b', ms=8)
#     h2, = ax.loglog(2./M, errs[1], 'v-g', ms=6)
#     h3, = ax.loglog(2./M, errs[2], 'x-r', ms=6)
#     h4, = ax.loglog(2./M, errs[3], 'd-', color='darkorange', ms=6)
#     h5, = ax.loglog(2./M, errs[4], 's-m', lw=2, ms=6)
#     h6, = ax.loglog(2./M, errsA[0], '.--b', ms=8)
#     h7, = ax.loglog(2./M, errsA[1], 'v--g', ms=6)
#     h8, = ax.loglog(2./M, errsA[2], 'x--r', ms=6)
#     h9, = ax.loglog(2./M, errsA[3], 'd--', color='darkorange', lw=2, ms=6)
#     h10, = ax.loglog(2./M, errsA[4], 's--m', lw=2, ms=6)
#     plt.xlabel('h',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Rel.error (L$_2$)', fontsize=18,fontweight='semibold')
#     ax.set(ylim=[1e-5,1e-1])
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':12}
#     plt.legend([h6,h7,h8,h9,h10,h1,h2,h3,h4,h5],[r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',], frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsBurgersGRKDG2hpCrConvAB'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# # ============================================
# with open(outputf+'BurgersGRKDGRKDG3hConvSL.npy', 'rb') as f:
#     errs = np.load(f)
# with open(outputf+'BurgersGRKDGRKDG2hConvSL.npy', 'rb') as f:
#     errsA = np.load(f)
# 
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.loglog(2./M, errs[0], '.-b', ms=8)
#     h2, = ax.loglog(2./M, errs[1], 'v-g', ms=6)
#     h3, = ax.loglog(2./M, errs[2], 'x-r', ms=6)
#     h4, = ax.loglog(2./M, errs[3], 'd-', color='darkorange', ms=6)
#     h5, = ax.loglog(2./M, errs[4], 's-m', lw=2, ms=6)
#     h6, = ax.loglog(2./M, errsA[0], '.--b', ms=8)
#     h7, = ax.loglog(2./M, errsA[1], 'v--g', ms=6)
#     h8, = ax.loglog(2./M, errsA[2], 'x--r', ms=6)
#     h9, = ax.loglog(2./M, errsA[3], 'd--', color='darkorange', lw=2, ms=6)
#     h10, = ax.loglog(2./M, errsA[4], 's--m', lw=2, ms=6)
#     plt.xlabel('h',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Rel.error (L$_2$)', fontsize=18,fontweight='semibold')
#     ax.set(ylim=[1e-5,1e-1])
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':12}
#     plt.legend([h6,h7,h8,h9,h10,h1,h2,h3,h4,h5],[r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',r'$\mathbf{\sigma = 0.1}$',r'$\mathbf{\sigma = 0.5}$',r'$\mathbf{\sigma = 1.0}$',\
#             r'$\mathbf{\sigma = 1.5}$',r'$\mathbf{\sigma = 2.0}$',], frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsBurgersGRKDG32hpCrConv'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
with open(outputf+'BurgersGRKDGRKDG3hConvSL.npy', 'rb') as f:
    errs = np.load(f)
    M = np.load(f)
    dt = np.load(f)

def createTable(errs, M, dt, sma):
    dterrs = errs*dt
    slopes = (np.log(dterrs[:,1:])-np.log(dterrs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
#     print('slopes =',slopes)
    ls = len(sma)
    v = ['' for k in range(ls+1)]
    print('\\begin{longtabu} to \\textwidth {r' + '|XX'.join(v) + '}')
    print('\t\\caption{Caption here...}')
    print('\t\\label{XXX:tab}\\tabularnewline')
    print('\t\\hline')
    print(('\t' + ' & \multicolumn{2}{c}{$\sigma = %1.1f$}'.join(v) + '\\centralize\\tabularnewline') % tuple(sma))
    print('\t\\hline')
    print('\t$M$' + ' & L$_2$ error'.join(v) + ' \\centralize \\tabularnewline')
    print('\t\\hline\\hline')
    w = [M[0]]
    for k in range(ls): w.append(errs[k,0])
    print(('\t$%3i$' + ' & $%5.5E$ & ---'.join(v) + ' \\tabularnewline')\
            % tuple(w))
    for k in range(len(M)-1):
        l = k+1
        w = [M[l]]
        for s in range(ls): w.append(errs[s,l]); w.append(slopes[s,k])
        print(('\t$%3i$' +' & $%5.5E$ & $%5.3f$'.join(v) + ' \\tabularnewline')\
                % tuple(w))
    print('\t\\hline')
    print('\\end{longtabu}')
    return None

sma = [0.1, 0.5, 1.0, 1.5, 2.0]
createTable(errs[:3], M, dt[:3], sma[:3])
createTable(errs[3:], M, dt[3:], sma[3:])
    
