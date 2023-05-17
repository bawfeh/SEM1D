# Numpy tutorial: basics
# import sys, os
import numpy as np
ext = '.eps'
# ext = '.pdf'
# folder = '/Users/bawfeh78/Documents/Research/Manuscripts2/Figures/'
folder = '/Users/bawfeh78/Documents/Research/Manuscripts3/Figures/'
outputf = './LinearRKDG/'

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, NullFormatter )
from matplotlib import ticker

# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# with open(outputf+'LinearRKDG1pConvCr2p5.npy', 'rb') as f:
#     errs = np.load(f)
#     qq = np.load(f)
# 
# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.semilogy(qq, errs[0], 'b.-', lw=2, ms=8)
#     h2, = ax.semilogy(qq, errs[1], 'rx-', lw=2, ms=6)
#     h3, = ax.semilogy(qq, errs[2], 'ms-', lw=2, ms=6)
#     plt.xlabel('q',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Error (L$_2$)', fontsize=18,fontweight='semibold')
#     # ax.set(xlim=[0,2.1])
#     ax.set(ylim=[1e-10,1e-0])
#     locmin = ticker.LogLocator(base=10, numticks=14)
#     ax.yaxis.set_minor_locator(locmin)
#     ax.yaxis.set_minor_formatter(ticker.NullFormatter())
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':14}
#     plt.legend([h1,h2,h3],['M = 20','M = 40', 'M = 80'], frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsLinearRKDG1hpConv'+ext, bbox_inches = 'tight')
#     plt.show()
# 
# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
with open(outputf+'LinearRKDG1hConvCr5.npy', 'rb') as f:
    errs = np.load(f)
    M = np.load(f)
    dt = np.load(f)

# for ext in ['.pdf','.eps']:
#     fig, ax = plt.subplots()
#     h1, = ax.loglog(2./M, errs[0], 'b.-', lw=2, ms=8)
#     h2, = ax.loglog(2./M, errs[1], 'rx-', lw=2, ms=6)
#     h3, = ax.loglog(2./M, errs[2], 'ms-', lw=2, ms=6)
#     plt.xlabel('h',fontsize=16,fontweight='semibold')
#     plt.ylabel(r'Error (L$_2$)', fontsize=18,fontweight='semibold')
#     ax.set(ylim=[1e-7,1e-0])
#     ax.tick_params(length=6, width=2, labelsize=14, direction='in')
#     ax.tick_params(which='minor', length=4, direction='in')
#     legend_props = {'weight':'bold','size':14}
#     plt.legend([h1,h2,h3],[r'$\mathbf{q = 1}$',r'$\mathbf{q = 2}$',r'$\mathbf{q = 3}$'],\
#             frameon=False,prop=legend_props,loc='best')
#     plt.savefig(folder+'errorsLinearRKDG1hpCr5Conv'+ext, bbox_inches = 'tight')
#     plt.show()
# # ===========================================

# errs *= dt
# slopes = (np.log(errs[:,1:])-np.log(errs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
# print('slopes =',slopes)


# ===========================================
def createTable(errs, M, dt):
    
    dterrs = errs*dt
    slopes = (np.log(dterrs[:,1:])-np.log(dterrs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
    print('slopes =',slopes)

    print('\\begin{longtabu} to \\textwidth {r|XX|XX|XX}')
    print('\t\\caption{Caption here...}')
    print('\t\\label{XXX:tab}\\tabularnewline')
    print('\t\\hline')
    print('\t  & \multicolumn{2}{|c|}{$q = 1$} & \multicolumn{2}{|c|}{$q = 2$} & \multicolumn{2}{|c|}{$q = 3$}\\centralize[3mm]\\tabularnewline')
    print('\t\\hline\\hline')
    print('\t$M$ & L$_2$ error & order & L$_2$ error & order  & L$_2$ error & order \\centralize[3mm] \\tabularnewline')
    print('\t\\hline')
    print('\t$%3i$ & $%5.5E$ & --- & $%5.5E$ & --- & $%5.5E$ & --- \\tabularnewline'\
            % (M[0], errs[0,0], errs[1,0], errs[2,0]))
    for k in range(len(errs[0])-1):
        l = k+1
        print('\t$%3i$ & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$ \\tabularnewline'\
                % (M[l], errs[0,l], slopes[0,k], errs[1,l], slopes[1,k], errs[2,l], slopes[2,k]))
    print('\t\\hline')
    print('\\end{longtabu}')
    return None

createTable(errs, M, dt)
    

# # ============================================
# fig, ax = plt.subplots()
# h1, = ax.loglog(2*np.pi/M, errs4[0], 'b.-', lw=2, ms=8)
# h2, = ax.loglog(2*np.pi/M, errs4[1], 'gv-', lw=2, ms=6)
# h3, = ax.loglog(2*np.pi/M, errs4[2], 'rx-', lw=2, ms=6)
# h4, = ax.loglog(2*np.pi/M, errs4[3], 'd-', lw=2, ms=6)
# h5, = ax.loglog(2*np.pi/M, errs4[4], 'ms-', lw=2, ms=6)
# plt.xlabel('h',fontsize=16,fontweight='semibold')
# plt.ylabel(r'Rel.error (L$_2$)', fontsize=18,fontweight='semibold')
# ax.set(ylim=[1e-7,1e+0])
# ax.tick_params(length=6, width=2, labelsize=14, direction='in')
# ax.tick_params(which='minor', length=4, direction='in')
# legend_props = {'weight':'bold','size':14}
# plt.legend([h1,h2,h3,h4,h5],[r'$\mathbf{\lambda = 0.0}$',r'$\mathbf{\lambda = 0.1}$',r'$\mathbf{\lambda = 0.5}$',\
#         r'$\mathbf{\lambda = 0.9}$',r'$\mathbf{\lambda = 1.0}$'], frameon=False,prop=legend_props,loc='lower right')
# plt.savefig('errorsLinearRKDG3hpCr3p5Convp3.'+ext, bbox_inches = 'tight')
# plt.show()
# 
# def createTable(M, errs):
#     slopes = (np.log(errs[:,1:])-np.log(errs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
#     
#     print('\\begin{longtabu} to \\textwidth {r|XX|XX|XX|XX|XX}')
#     print('\t\\caption{Caption here...}')
#     print('\t\\label{XXX:tab}\\tabularnewline')
#     print('\t\\hline')
#     print('\t  & \multicolumn{2}{c}{$\lambda = %1.1f$} & \multicolumn{2}{|c}{$\lambda = %1.1f$} & \multicolumn{2}{|c}{$\lambda = %1.1f$} & \multicolumn{2}{|c}{$\lambda = %1.1f$} & \multicolumn{2}{|c}{$\lambda = %1.1f$}\\centralize\\tabularnewline' % (0.0,0.1,0.5,0.9,1.0))
#     print('\t\\hline')
#     print('\t$M$ & L$_2$ error & order & L$_2$ error & order  & L$_2$ error & order',\
#             '& L$_2$ error & order  & L$_2$ error & order \\centralize \\tabularnewline')
#     print('\t\\hline\\hline')
#     print('\t$%3i$ & $%5.5E$ & --- & $%5.5E$ & --- & $%5.5E$ & --- & $%5.5E$ & --- & $%5.5E$ & --- \\tabularnewline'\
#             % (M[0], errs[0][0], errs[1][0], errs[2][0], errs[3][0], errs[4][0]))
#     for k in range(len(M)-1):
#         l = k+1
#         print('\t$%3i$ & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$ \\tabularnewline'\
#                 % (M[l], errs[0][l], slopes[0][k], errs[1][l], slopes[1][k], errs[2][l], slopes[2][k], errs[3][l], slopes[3][k], errs[4][l], slopes[4][k]))
#     print('\t\\hline')
#     print('\\end{longtabu}')
#     return None
# 
# def createTable(M, errs):
#     slopes = (np.log(errs[:,1:])-np.log(errs[:,:-1])) / (np.log(M[:-1])-np.log(M[1:]))
#     
#     print('\t\\hline\\hline')
#     print('\t$%3i$ & $%5.5E$ & --- & $%5.5E$ & --- & $%5.5E$ & --- \\tabularnewline'\
#             % (M[0], errs[0][0], errs[1][0], errs[2][0]))
#     for k in range(len(M)-1):
#         l = k+1
#         print('\t$%3i$ & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  & $%5.5E$ & $%5.3f$  \\tabularnewline'\
#                 % (M[l], errs[0][l], slopes[0][k], errs[1][l], slopes[1][k], errs[2][l], slopes[2][k]))
#     print('\t\\hline')
#     print('\t\\hline\\hline')
#     print('\t$%3i$ & $%5.5E$ & --- & $%5.5E$ & \multicolumn{1}{l}{ --- } \\tabularnewline'\
#             % (M[0], errs[3][0], errs[4][0]))
#     for k in range(len(M)-1):
#         l = k+1
#         print('\t$%3i$ & $%5.5E$ & $%5.3f$  & $%5.5E$ & \multicolumn{1}{l}{ $%5.3f$ } \\tabularnewline'\
#                 % (M[l], errs[3][l], slopes[3][k], errs[4][l], slopes[4][k]))
#     print('\t\\hline')
#     print('\\end{longtabu}')
#     return None
# 
# createTable(M, errs4)
#   
