import matplotlib.pyplot as plt
import numpy as np
import plot_lib as pl
import matplotlib.cm as cm

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

def calc_heatmap(filename_id, energy):
    f = open(filename_id + '_' + str(energy) + '_heatmap.txt', 'wb')

    tt = np.linspace(0, 300, 151)
    rr = np.linspace(0, 200, 151)
    
    for t_ in tt:
        filename = 'output/fcr_' + filename_id + '_t_' + str(int(t_)) + '_nz_1201_np_256.txt'
        print filename
        x, y = read_spectrum_at_pz(filename, 3, energy)

        for r_ in rr:
            idx = (np.abs(x - r_)).argmin()
            #dzz.append(np.log10(y[idx] / D_ISM))
            f.write("%.4e   %.4e   %.4e\n" % (t_, r_, y[idx]))

#calc_heatmap('new_geminga_0', 1e2)
#calc_heatmap('new_geminga_0', 1e4)

#calc_heatmap('new_geminga_5', 1e2)
#calc_heatmap('new_geminga_5', 1e4)

#calc_heatmap('new_geminga_6', 1e2)
#calc_heatmap('new_geminga_6', 1e4)

#calc_heatmap('new_geminga_7', 1e2)
#calc_heatmap('new_geminga_7', 1e4)

#calc_heatmap('new_geminga_1', 1e2)
#calc_heatmap('new_geminga_1', 1e4)

#calc_heatmap('new_geminga_2', 1e2)
#calc_heatmap('new_geminga_2', 1e4)

#calc_heatmap('new_geminga_3', 1e2)
#calc_heatmap('new_geminga_3', 1e4)

#############################################
fig, ax = pl.set_plot_style()

filename = 'new_geminga_0_10000.0_heatmap.txt'
plot_filename = 'D_0_1e4_GeV.png'
energy = 1e4
#############################################

D_ISM = 5e28 * (energy / 3.)**(1./3.)

ax.set_ylabel(r't [kyr]', fontsize = 30, labelpad = 10)
ax.set_xlabel(r'r [pc]', fontsize = 30, labelpad = 15)
#ax.set_title(plot_title, fontsize = 30)

t_, r_, dzz_ = np.loadtxt(filename, skiprows=0, usecols=(0,1,2), unpack=True)

print len(dzz_)

tt = t_[0:151*151:151]
rr = r_[0:151]

log_dzz = np.reshape(np.log10(dzz_ / D_ISM), [len(rr), len(tt)])

##im = ax.pcolormesh(rr, tt, dzz, vmin=27, vmax=30, cmap='jet_r')

##cmap='RdBu' coolwarm RdYlBu
#clb = fig.colorbar(im)
#clb.set_label(r'log D / D$_{ISM}$', fontsize = 28)

im = plt.imshow(log_dzz, interpolation='bilinear', cmap=cm.winter, origin='lower', extent=[0, 100, 0, 200], vmax=0, vmin=-3)

CS = plt.contour(tt, rr, log_dzz, levels=[-2, -1.3, -1, -0.3, 0], colors='white', linestyles='-', linewidhts=8)

fmt = {}
strs = [r'1\%', r'5\%', r'10\%', '50\%', '1']
for l, s in zip(CS.levels, strs):
    fmt[l] = s

plt.clabel(CS, CS.levels[::1], fmt=fmt, fontsize=25, inline=1) 
#plt.clabel(CS, inline=1, fontsize=25, fmt=fmt)

#im = ax.pcolormesh(rr, tt, log_dzz, vmin=-3, vmax=0, cmap='jet_r')

#plt.text(60, 90, 'E = 10 TeV', fontsize=27, color='w')

ax.set_xlim([0,100])
ax.axis('tight')

plt.xlim([0, 100])

#plt.show()

plt.savefig(plot_filename)

#r_positron = 0.0
#t_positron = 40
#dt_positron = 1e3 * 3.14e7 # 1 kyr
#if t_ > t_positron:
#        idx = (np.abs(x - r_positron)).argmin()
#       dzz_positron = y[idx]
#        delta_r = np.sqrt(dzz_positron * dt_positron) / 3.1e18
#        r_positron += delta_r
#        if r_positron < 100:
#            r_p.append(r_positron)
#            t_p.append(t_)
#            print t_positron / 3.14e7, dzz_positron, delta_r, r_positron
#ax.plot(r_p, t_p, 'white')
#ax.set_xlim([0,100])
#ax.set_ylim([0,100])

