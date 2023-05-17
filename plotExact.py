import pickle
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

ext = 'eps'
source = '/Users/bawfeh78/Documents/PYTHON/SEM1D/Data/'
dest = ''

with open(source+'odeRadau.pickle', 'rb') as f:
    exactsol = pickle.load(f)

from matplotlib import pyplot as plt
fig, ax = plt.subplots()
h1, = ax.plot(exactsol[-1].nodes, exactsol[0][:,0], linewidth=2)
h2, = ax.plot(exactsol[-1].nodes, exactsol[0][:,81], '--', linewidth=2)
h3, = ax.plot(exactsol[-1].nodes, exactsol[0][:,120], ':', linewidth=2)
h4, = ax.plot(exactsol[-1].nodes, exactsol[0][:,139], '-.', linewidth=2)
h5, = ax.plot(exactsol[-1].nodes, exactsol[0][:,-1], '.-', linewidth=2, color='purple', markevery=80)
plt.xlabel('x',fontsize=18,fontweight='bold')
# plt.axis([-.05,1.01,-.01,1.02])
ax.tick_params(length=6, width=2, labelsize=14,direction='in')
ax.set(xlim=[-0.1,1.01])
ax.set(ylim=[0.0,1.01])
# ax.xaxis.grid(True, which='minor')
#     ax.yaxis.set_major_locator(MultipleLocator(.1))
#     ax.yaxis.set_minor_locator(MultipleLocator(.05))
#     ax.xaxis.set_major_locator(MultipleLocator(.5))
#     ax.xaxis.set_minor_locator(MultipleLocator(.1))
legend_props = {'weight':'bold','size':14}
plt.legend([h1,h2,h3,h4,h5],['t=0.0','t=0.5026','t=1.0016','t=1.5177','t=2.0'], frameon=False,prop=legend_props,loc='best')
plt.savefig(dest+'exactBurgers1D.'+ext, bbox_inches = 'tight')#, pad_inches = 0)
plt.show()

# ============================================
