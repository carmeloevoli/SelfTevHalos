import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def read_profile_at_p(filename, column, p_search, linestyle, donormalize):
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
    print "plot f in the range = ", min(f), max(f), np.mean(f)
    if donormalize:
        plt.plot(z, f / max(f), linestyle)
    else:
        plt.plot(z, f, linestyle)

def read_spectrum_at_z(filename, column, z_search, pslope, linestyle):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    f = f[ind] # GeV^-1 m^-2 s^-1
    print "plot f in the range = ", min(f), max(f), np.mean(f)
    plt.plot(p, p**pslope * f, linestyle)

def get_spectrum_at_z(filename, column, z_search):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    f = f[ind] # GeV^-1 m^-2 s^-1
    return p, f
    
def plot_data(filename, slope_plot, color, label):
    T, y, err_y_lo, err_y_up = np.loadtxt(filename, skiprows=2, usecols=(0,3,8,9), unpack=True)
    y_err = [T**slope_plot * err_y_lo, T**slope_plot * err_y_up]
    plt.errorbar(T, T**slope_plot * y, yerr=y_err, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color=color)
    
def set_plot_style():
    #plt.style.use('bmh')
    #print(matplotlib.rcParams.keys())
    matplotlib.rcParams.update({
        #'axes.grid': True,
        #'axes.titlesize': 'medium',
        'font.family': 'serif',
        'font.serif': 'Helvetica Neue',
        'font.size': 30,
        #'grid.color': 'w',
        #'grid.linestyle': '-',
	#'grid.alpha': 0.5,
 	#'grid.linewidth': 1,
 	'legend.frameon': False,
 	'legend.fancybox': False,
 	'legend.fontsize': 20,
 	#'legend.framealpha': 0.7,
 	#'legend.handletextpad': 0.1,
 	#'legend.labelspacing': 0.2,
 	'legend.loc': 'best',
 	'lines.linewidth': 3,
 	'savefig.bbox': 'tight',
 	#'savefig.pad_inches': 0.02,
        'figure.autolayout': True,
 	'text.usetex': True, 
 	#'text.latex.preamble': r'\usepackage{txfonts}',
        'xtick.labelsize': 30,
        'ytick.labelsize': 30,
        'axes.labelpad': 10,
 	})
    fig = plt.figure(figsize=(9.0, 8.6))
    ax = fig.add_subplot(1, 1, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.minorticks_on()
    ax.tick_params('both', length=15, width=1.5, which='major', pad=10)
    ax.tick_params('both', length=0,  width=1.3, which='minor', pad=10)
    return fig, ax

matplotlib.rc("savefig", dpi=200)
matplotlib.rc("figure", figsize=(9.0, 8.67))


# PLOT STYLE

#rcParams['legend.numpoints'] = 1
#rcParams['lines.linewidth'] = 3
#rcParams['figure.autolayout'] = True

#fig = plt.figure(figsize=(9.0, 8.6))
#ax = fig.add_subplot(1, 1, 1)

#plt.xticks(size=35)
#plt.yticks(size=35)
# END PLOT STYLE
