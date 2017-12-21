#!/bin/bash/python
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np

# PLOT STYLE
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='Helvetica Neue')
rc('xtick', labelsize=18)
rc('ytick', labelsize=18)
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 4
rcParams['figure.autolayout'] = True

fig = plt.figure(figsize=(9.2, 8.6))
ax = fig.add_subplot(1, 1, 1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.minorticks_on()
ax.tick_params('both', length=15, width=1.5, which='major', pad=6)
ax.tick_params('both', length=10, width=1.3, which='minor', pad=6)

plt.xticks(size=30)
plt.yticks(size=30)
# END PLOT STYLE

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def read_spectrum_at_pz(filename, column, p_search):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (p_search < min(p) or p_search > max(p)):
        print "WARNING: p is out of range!"
    best_p = -1000.0
    for i in range(len(p)):
        if abs(p[i] - p_search) < abs(best_p - p_search):
            best_p = p[i]
    ind = (p == best_p).nonzero()
    z = z[ind]
    f = f[ind]
    return z, f

energy = 1e2
plot_title = '100 GeV'
plot_filename = 'heatmap_100_GeV.png'
D_ISM = 5e28 * (energy / 3.)**(1./3.)

ax.set_ylabel(r't [kyr]', fontsize = 30, labelpad = 10)
ax.set_xlabel(r'r [pc]', fontsize = 30, labelpad = 15)
ax.set_title(plot_title, fontsize = 30)

tt = np.linspace(0, 340, 171)
rr = np.linspace(0,  50, 171)

r_positron = 0.0
t_positron = 40
dt_positron = 1e3 * 3.14e7 # 1 kyr

dzz = []
r_p = []
t_p = []

for t_ in tt:
    filename = 'output/fcr_test_phase_8_t_' + str(int(t_)) + '_nz_1001_np_160.txt'
    print filename
    x, y = read_spectrum_at_pz(filename, 3, energy)

    for r_ in rr:
        idx = (np.abs(x - r_)).argmin()
        #print r_, x[idx], y[idx]
        dzz.append(np.log10(y[idx]))

#if t_ > t_positron:
#        idx = (np.abs(x - r_positron)).argmin()
#       dzz_positron = y[idx]
#        delta_r = np.sqrt(dzz_positron * dt_positron) / 3.1e18
#        r_positron += delta_r
#        if r_positron < 100:
#            r_p.append(r_positron)
#            t_p.append(t_)
#            print t_positron / 3.14e7, dzz_positron, delta_r, r_positron

dzz = np.reshape(dzz, [len(rr), len(tt)])

#fig, ax = plt.subplots()
im = ax.pcolormesh(rr, tt, dzz, vmin=25.5, vmax=29.5, cmap='jet_r')
#cmap='RdBu' coolwarm RdYlBu
clb = fig.colorbar(im)
clb.set_label(r'log D [cm$^2$ s$^{-1}$]', fontsize = 28)
#ax.plot(r_p, t_p, 'white')
#ax.set_xlim([0,100])
#ax.set_ylim([0,100])

ax.axis('tight')

#plt.show()
plt.savefig(plot_filename, format='png', dpi=300)
