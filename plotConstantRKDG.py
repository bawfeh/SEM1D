# Numpy tutorial: basics
# import sys, os
import numpy as np
# ext = 'eps'
# ext = 'pdf'
# folder = '/Users/bawfeh78/Documents/Research/Manuscripts2/Figures/'
folder = '/Users/bawfeh78/Documents/Research/Manuscripts3/Figures/'
outputf = './ConstantRKDG/'

from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator, NullFormatter )
from matplotlib import ticker

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

sol = []; err = []
with open(outputf+'ShuSLRKDG1A.npy', 'rb') as f:
    sol.append(np.load(f))
    x = np.load(f)
    t = np.load(f)
    exactsol = np.load(f)
    err.append(np.load(f))

with open(outputf+'ShuSLRKDG1B.npy', 'rb') as f:
    sol.append(np.load(f))
    x = np.load(f)
    t = np.load(f)
    exactsol = np.load(f)
    err.append(np.load(f))

# print('x =', x)
print('Error = ', err)
print('t =', t)

for ext in ['.pdf','.eps']:
    figname = [folder+'ShuSLRKDG1B'+ext,folder+'ShuSLRKDG1B2'+ext]
    
    for k in range(2):
        fig, ax = plt.subplots()
        h1, = ax.plot(x, exactsol, lw=2, ms=6)
        h2, = ax.plot(x, sol[k],'r.--', lw=1, ms=6)
        ax.set(ylim=[-.05,1.05])
        plt.xlabel('x',fontsize=16,fontweight='semibold')
        plt.ylabel('u',fontsize=16,fontweight='semibold')
        ax.yaxis.set_major_locator(MultipleLocator(.2))
        ax.yaxis.set_minor_locator(MultipleLocator(.1))
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.xaxis.set_minor_locator(MultipleLocator(.1))
        ax.tick_params(length=6, width=2, labelsize=14, direction='in')
        ax.tick_params(which='minor', length=4, direction='in')
        legend_props = {'weight':'bold','size':14}
        plt.legend([h1,h2],['Exact','SLDG'], frameon=False,prop=legend_props,loc='upper right')
        print(figname[k])
        plt.savefig(figname[k], bbox_inches = 'tight')
        plt.show()

# print('\\begin{figure}[t]%[thbp]')
# print('  \centering')
# print('  \subfigure[]{')
# print('    \includegraphics[width=0.45\\textwidth]{ShuSLRKDG1B.pdf}')
# print('    \label{fig:ShuSLRKDG1B}')
# print('  }%~\hspace{-1.0cm}')
# print('  \subfigure[]{')
# print('    \includegraphics[width=0.45\\textwidth]{ShuSLRKDG1B2.pdf}')
# print('    \label{fig:ShuSLRKDG1B2}')
# print('  }')
# print('  \caption{}')
# print('  \label{fig:ShuSLRKDG}')
# print('\end{figure}')
# 
