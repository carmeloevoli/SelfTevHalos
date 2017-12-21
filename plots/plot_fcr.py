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
rcParams['lines.linewidth'] = 3
rcParams['figure.autolayout'] = True

fig = plt.figure(figsize=(9.0, 8.6))
ax = fig.add_subplot(1, 1, 1)

for axis in ['top', 'bottom', 'left', 'right']:
    ax.spines[axis].set_linewidth(1.5)

ax.minorticks_on()
ax.tick_params('both', length=15, width=1.5, which='major', pad=10)
ax.tick_params('both', length=0,  width=1.3, which='minor', pad=10)

plt.xticks(size=35)
plt.yticks(size=35)
# END PLOT STYLE

def plot_data(filename, slope_plot, color, label, norm):
    T, y, err_y_lo, err_y_up = np.loadtxt(filename,skiprows=2,usecols=(0,3,8,9),unpack=True)
    y_err = [T**slope_plot * err_y_lo, T**slope_plot * err_y_up]
    plt.errorbar(T, T**slope_plot * y, yerr=y_err, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color=color)

def read_spectrum_at_pz(filename, column, z_search, p_search):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (p_search < min(p) or p_search > max(p)):
        print "WARNING: p is out of range!"
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: p is out of range!"
    best_p = -1000.0
    for i in range(len(p)):
        if abs(p[i] - p_search) < abs(best_p - p_search):
            best_p = p[i]
    ind = (p == best_p).nonzero()
    z = z[ind]
    f = f[ind]

    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    f = f[ind]
    print 'f is',f
    return f

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

def read_spectrum_at_z(filename, column, z_search, alpha, linestyle):
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
    plt.plot(p, p**alpha * f, linestyle)

def read_energy_density_at_z(filename, column, z_search):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    f = f[ind] # GeV/c^-3 m^-3

    epsilon = np.sum(p**4 * f)
    dlnp = np.log(p[1] / p[0])
    
    return 4.0 * 3.14 * epsilon * dlnp * 1e3



def read_v_at_z(filename, z_search, alpha, linestyle):
    z, p, f, D, dfdz = np.loadtxt(filename, skiprows=1, usecols=(0,1,2,3,4), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    f = f[ind] # GeV^-1 m^-2 s^-1
    D = D[ind]
    dfdz = dfdz[ind]
    plt.plot(p, p**alpha * D / f * dfdz / 3e21, linestyle)

def read_timescale(filename, n_z, n_p, alpha, color):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,3), unpack=True)
    dz = (max(z) - min(z)) / n_z
    print max(z), dz
    counter = 0
    p_x = np.zeros(n_p)
    I_y = np.zeros(n_p)
    I_small_y = np.zeros(n_p)
    for iz in range(n_z):
        for ip in range(n_p):
            p_x[ip] = p[counter]
            I_y[ip] += abs(z[counter]) / f[counter]
            if abs(z[counter]) > 2.0:
                I_small_y[ip] += abs(z[counter]) / f[counter]
            #print counter, z[counter]
            counter = counter + 1
    plt.plot(p_x, p_x**alpha * dz * I_y * 31., color=color)
    plt.plot(p_x, p_x**alpha * dz * I_small_y * 31., color=color, linestyle=':')

def read_jcr_at_z(filename, column, z_search, alpha, linestyle):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    E = p # GeV
    f = f[ind] # m^-2 s^-1
    EJ = f
    print "plot J in the range = ", min(EJ), max(EJ),np.mean(EJ)
    plt.plot(E, E**alpha * EJ, linestyle)

def analytical_solution_profile(H, D, vA): # kpc - 1e28 cm2/s - km/s
    xi = 0.03 * vA * H / D
    z = np.linspace(-H,H,100)
    f = (1. - np.exp(-xi * (1. - abs(z) / H))) / (1. - np.exp(-xi))
    plt.plot(z, f, 'k:')

def kolmogorov_coefficient(D0, delta):
    m_p = 0.938
    p = np.logspace(-2,10,120)
    E = np.sqrt(p * p + m_p * m_p)
    beta = p / E
    plt.plot(p, D0 * beta * p**delta, 'k:')

def advection_timescale(vA, H):
    p = np.logspace(-3,10,100)
    t = ((H * 3.086e+19) / (vA * 1e3))  * (p / p)
    plt.plot(p, t / 3.14e13, 'k:')

def read_ratio_at_z(fn0, fn, column, z_search, alpha, linestyle):
    z, p, f0 = np.loadtxt(fn0, skiprows=1, usecols=(0,1,column), unpack=True)
    z, p, f  = np.loadtxt(fn,  skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print "WARNING: z is out of range!"
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p  = p[ind] # GeV/c
    f0 = f0[ind] # GeV^-1 m^-2 s^-1
    f  = f[ind]
    print "plot f in the range = ", min(f), max(f)
    print "plot f0 in the range = ", min(f0), max(f0)
    plt.plot(p, (f - f0) / f0, linestyle)

def plot_Jcr_spectrum():
    alpha = 2
    #plot_data('data/H_BESS-PolarII.txt', alpha, 'c', 'BESS-PolarII', 1)
    #plot_data('data/H_PAMELA.txt', alpha, '#FFA500', 'PAMELA', 1)
    #plot_data('data/H_AMS02.txt', alpha, 'r', 'AMS-02', 1)
    #plot_data('data/H_CREAM.txt', alpha, 'g', 'CREAM', 1)
    #plot_data('data/H_Voyager.txt', alpha, 'm', 'Voyager', 1)

#plt.legend(['BESS-PolarII', 'PAMELA', 'AMS-02', 'CREAM', 'Voyager'], fontsize=20, loc='upper left', frameon=False)

    plt.xscale('log')
    plt.xlim([1e0,1e5])

    plt.xlabel(r'T [GeV]', fontsize=38)
    
    plt.yscale('log')
    #plt.ylim([1e2, 1e5])
    plt.ylabel(r'T$^{2}$ J(T) [GeV/(m$^2$ s sr)]', fontsize=38)
    
    filename = 'output/fcr_tevhalo_3_t_100_nz_401_np_160.txt'
    read_jcr_at_z(filename, 2, 0.1, alpha, 'orange', 1)
    read_jcr_at_z(filename, 2, 0.2, alpha, 'g', 1)
    read_jcr_at_z(filename, 2, 0.4, alpha, 'r', 1)

    return 'electron_spectrum.pdf'

def plot_fcr_spectrum():
    alpha = 3
    p = np.logspace(-2, 6, 100)
    
    filename = 'output/fcr_test_1D_1_t_1_nz_1001_np_20.txt'
    read_spectrum_at_z(filename, 2,    5, alpha, 'y')
    read_spectrum_at_z(filename, 2,   10, alpha, 'g')
    read_spectrum_at_z(filename, 2,   20, alpha, 'r')

    plt.xscale('log')
    #plt.xlim([1e2, 1e5])
    plt.xlabel(r'p [GeV/c]',fontsize=30)
    
    plt.yscale('log')
    plt.ylim([1e-15, 1e-5])
    plt.ylabel(r'p$^3$ f(p) []',fontsize=30)

    p = np.logspace(2, 5, 100)
    #plt.plot(p, 0.5e-11 * (p / 1e4)**(-3.5), 'k:')
    #plt.plot(p, 1.0e-7 * (p / 1e2)**(-1.5), 'k:')
    return 'fcr.pdf'

def plot_fcr_profile():
    
    filename = 'output/fcr_test_3D_0_t_1_nz_401_np_20.txt'
    read_profile_at_p(filename, 2, 1e5, 'r', True)
    read_profile_at_p(filename, 2, 1e4, 'b', True)
    read_profile_at_p(filename, 2, 1e3, 'm', True)

    filename = 'output/fcr_test_3D_0_t_2_nz_401_np_20.txt'
    read_profile_at_p(filename, 2, 1e5, 'r:', True)
    read_profile_at_p(filename, 2, 1e4, 'b:', True)
    read_profile_at_p(filename, 2, 1e3, 'm:', True)
    
#    filename = 'output/fcr_test_1D_1_t_4_nz_1001_np_20.txt'
#    read_profile_at_p(filename, 2, 1e5, 'r:', True)
#    read_profile_at_p(filename, 2, 1e3, 'b:', True)
    
    #plt.xscale('log')
    plt.xlim([0, 10])
    plt.xlabel(r'z [pc]', fontsize=28)
    
    #plt.yscale('log')
    #plt.ylim([0.8, 1])
    plt.ylabel(r'f(z) [au]', fontsize=28)
    
    #z = np.logspace(-1, 3, 100)
    #plt.plot(z, np.exp(-0.5 * (z / 25.)**2), 'k--')

#plt.legend(['10 GeV', '100 GeV', '1 TeV', '10 TeV'], fontsize=25, loc='upper right', frameon=False)

    return 'fcr_profile.pdf'

def plot_giovanni():

    filename = 'output/fcr_geminga_0_t_300_nz_3001_np_160.txt'
    read_v_at_z(filename, 5,  0, 'y')
    read_v_at_z(filename, 10, 0, 'g')
    read_v_at_z(filename, 20, 0, 'r')
    plt.xscale('log')
    plt.yscale('log')


def plot_dzz_profile():

    filename = 'output/fcr_test_1D_1_t_59_nz_1001_np_20.txt'
    read_profile_at_p(filename, 3, 1e2, 'y', False)
    read_profile_at_p(filename, 3, 1e3, 'g', False)
    read_profile_at_p(filename, 3, 1e4, 'r', False)
    read_profile_at_p(filename, 3, 1e5, 'b', False)

    plt.xlabel(r'z [kpc]', fontsize=38)
    #plt.xlim([0,10])
    #plt.xscale('log')

    plt.ylabel(r'D(z) [$10^{28}$ cm$^2$/s]', fontsize=38)
    #plt.ylim([0.01,1.1])
    plt.yscale('log')

    plt.legend([r'$E = 100$ TeV'], fontsize=25, loc='upper right', frameon=False)

    return 'dzz_profile.pdf'

def plot_dzz_withtime():
    t = []
    d0 = []
    d1 = []
    d2 = []
    for i in range(105):
        t.append(i)
        filename = 'output/fcr_geminga_0_t_' + str(i) + '_nz_3001_np_160.txt'
        dzz = read_spectrum_at_pz(filename, 3, 0.01, 1e4)
        d0.append(dzz)
        dzz = read_spectrum_at_pz(filename, 3, 0.05, 1e4)
        d1.append(dzz)
        dzz = read_spectrum_at_pz(filename, 3, 0.25, 1e4)
        d2.append(dzz)
    
    plt.plot(t, d0, 'g')
    plt.plot(t, d1, 'r')
    plt.plot(t, d2, 'b')

    plt.xlim([0,100])
    #plt.ylim([1e26, 1e29])

    plt.yscale('log')

    plt.legend(['10 pc', '50 pc', '250 pc'], fontsize=25, loc='upper right', frameon=False)

    return 'dzz_with_time.pdf'

def plot_dzz_spectrum():
    alpha = 0 # -1. / 3.

    filename = 'output/fcr_test_1D_1_t_2_nz_101_np_20.txt'
    #read_spectrum_at_z(filename, 3,    -1, alpha, 'r:')
    #read_spectrum_at_z(filename, 3,  -0.2, alpha, 'b:')
    #read_spectrum_at_z(filename, 3,  -0.1, alpha, 'g:')
    #read_spectrum_at_z(filename, 3,     0, alpha, 'y')
    #read_spectrum_at_z(filename, 3,   0.1, alpha, 'g')
    #read_spectrum_at_z(filename, 3,   0.2, alpha, 'b')
    read_spectrum_at_z(filename, 3,     1, alpha, 'r')
    
    filename = 'output/fcr_test_1D_1_t_2_nz_201_np_20.txt'
    read_spectrum_at_z(filename, 3,     1, alpha, 'b')

    filename = 'output/fcr_test_1D_1_t_2_nz_401_np_20.txt'
    read_spectrum_at_z(filename, 3,     1, alpha, 'm')

    filename = 'output/fcr_test_1D_1_t_2_nz_801_np_20.txt'
    read_spectrum_at_z(filename, 3,     1, alpha, 'c')
    
    #read_spectrum_at_z(filename, 3,    10, alpha, 'm')
    #read_spectrum_at_z(filename, 3,    50, alpha, 'm')

    #filename = 'output/fcr_test_1D_1_t_0_nz_2001_np_20.txt'
    #read_spectrum_at_z(filename, 3,  1, alpha, 'g:')
    #read_spectrum_at_z(filename, 3, 10, alpha, 'r:')
    #read_spectrum_at_z(filename, 3, 20, alpha, 'b:')

    #plt.legend(['10 pc', '20 pc', '50 pc'], fontsize=22, loc='upper left', frameon=False)

    #filename = 'output/fcr_test_1D_1_t_0_nz_2001_np_20.txt'
    #read_spectrum_at_z(filename, 3, 10, alpha, 'k:')
    
    x = 2e4
    y = 3.2e27 * (100. / 20.)**(1./3.)
    y_err_hi = 1.4e27 * (100. / 20.)**(1./3.)
    y_err_lo = 1.0e27 * (100. / 20.)**(1./3.)

    #plt.errorbar(x, x**alpha * y, yerr=y_err_lo, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color='k')

    #plt.plot([5e4, 5e4], [1e20, 1e40], 'c:')

    p = np.logspace(2,5,100)
    
    c = 3e10 # cm / s
    rL = 3.1e18 * (p / 1e6)
    
    plt.plot(p, rL * c / 3.0, 'k--')
    
    #plt.text(1e4, 2e27, 'H', fontsize=20)
    
    plt.text(1e3, 2e30, 'Galactic Diffusion', fontsize=25, rotation=11)
    plt.text(1e4, 5e26, 'Bohm limit', fontsize=25, rotation=30)

    plt.xscale('log')
    #plt.xlim([1e1, 1e6])
    plt.xlabel(r'p [GeV/c]',fontsize=35,labelpad=10)

    plt.yscale('log')
    #plt.ylim([1e26, 1e31])
    plt.ylabel(r'D(p) [cm$^2$/s]',fontsize=35,labelpad=15)

    return 'dzz_spectrum.pdf'

def plot_dzz_timescale():
    alpha = 0

    read_spectrum_at_z('output/fcr_PRL_1_t_0_nz_257_np_704.txt', 4, 0.0, alpha, 'g')
    read_spectrum_at_z('output/fcr_PRL_1_t_0_nz_257_np_704.txt', 4, 0.2, alpha, '#FFA500')
    read_spectrum_at_z('output/fcr_PRL_1_t_0_nz_257_np_704.txt', 4, 0.5, alpha, 'r')
    read_spectrum_at_z('output/fcr_PRL_1_t_0_nz_257_np_704.txt', 4, 1.0, alpha, 'b')
    read_spectrum_at_z('output/fcr_PRL_1_t_0_nz_257_np_704.txt', 4, 2.0, alpha, 'm')

    #read_spectrum_at_z('output/fcr_BAS_0_t_0_nz_257_np_600.txt', 4, 0.0, alpha, 'r')
    #read_spectrum_at_z('output/fcr_BAS_0_t_1_nz_321_np_320.txt', 4, 0.0, alpha, 'b:')
    #read_spectrum_at_z('output/fcr_BAS_0_t_2_nz_321_np_320.txt', 4, 0.0, alpha, 'g--')
    
    #read_spectrum_at_z('output/fcr_BAS_0_dt_7.55948e-06_t_0_nz_321_np_320.txt', 4, 0, alpha, 'r')
    #read_spectrum_at_z('output/fcr_BAS_1_dt_7.55948e-06_t_0_nz_321_np_320.txt', 4, 0, alpha, 'b')
    #read_spectrum_at_z('output/fcr_BAS_0_dt_7.55948e-06_t_0_nz_321_np_288.txt', 4, 0.1, alpha, 'r')
    #read_spectrum_at_z('output/fcr_BAS_0_dt_7.55948e-06_t_0_nz_321_np_288.txt', 4, 1.0, alpha, 'r')
    #read_spectrum_at_z('output/fcr_BAS_0_dt_7.55948e-06_t_0_nz_321_np_288.txt', 4, 2.0, alpha, 'r')
    advection_timescale(10.0, 4.0)

    plt.xscale('log')
    plt.xlim([0.01,1e8])
    plt.xlabel(r'p [GeV/c]',fontsize=28)

    plt.yscale('log')
    plt.ylim([1e-2,1e5])
    plt.ylabel(r't$_{zz}$(p) [Myr]',fontsize=28)

def plot_dfdz_spectrum():
    alpha = 3.0

    filename = 'output/fcr_test_1D_1_t_1_nz_1001_np_20.txt'
    read_spectrum_at_z(filename, 4, -1, alpha, 'g')
    read_spectrum_at_z(filename, 4, -2, alpha, 'c')
    read_spectrum_at_z(filename, 4, -5, alpha, 'b')

    read_spectrum_at_z(filename, 4,  0, alpha, 'y')
    
    read_spectrum_at_z(filename, 4,  1, alpha, 'g:')
    read_spectrum_at_z(filename, 4,  2, alpha, 'c:')
    read_spectrum_at_z(filename, 4,  5, alpha, 'b:')

    plt.xscale('log')
    #plt.xlim([0.1,1e7])
    plt.xlabel(r'p [GeV/c]',fontsize=28)
    
    plt.yscale('log')
    #plt.ylim([1e-20,1e0])
    plt.ylabel(r'p3 df/dz []',fontsize=28)

def plot_dfdz_profile():
    filename = 'output/fcr_test_1D_1_t_3_nz_1001_np_20.txt'
    read_profile_at_p(filename, 4, 1e2, 'g', True)
    read_profile_at_p(filename, 4, 1e3, 'c', True)
    read_profile_at_p(filename, 4, 1e4, 'b', True)
    read_profile_at_p(filename, 4, 1e5, 'm', True)

    #plt.xscale('log')
    plt.xlim([-5, 5])
    plt.xlabel(r'z [kpc]',fontsize=28)
    
    plt.yscale('log')
    #plt.ylim([1e-5,1e1])
    plt.ylabel(r'p$^3$ df/dz []',fontsize=28)

    return 'dfdz_profile.pdf'

def plot_flux_profile():
    filename = 'output/fcr_test_1D_0_t_10_nz_1001_np_20.txt'
    read_profile_at_p(filename, 5, 1e2, 'g', False)
    read_profile_at_p(filename, 5, 1e3, 'c', False)
    read_profile_at_p(filename, 5, 1e4, 'b', False)
    read_profile_at_p(filename, 5, 1e5, 'm', False)
    
    plt.xlim([90, 100])
    plt.yscale('log')

def plot_wave_profile():
    filename = 'output/fcr_test_1D_1_t_10_nz_101_np_20.txt'
    read_profile_at_p(filename, 6, 1e2, 'y', True)
    read_profile_at_p(filename, 6, 1e3, 'g', True)
    read_profile_at_p(filename, 6, 1e4, 'r', True)
    read_profile_at_p(filename, 6, 1e5, 'b', True)

    plt.xlim([-10, 10])
    return 'wave_profile.pdf'


def plot_Gamma_spectrum():
    alpha = 0
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.0, alpha, 'y')
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.1, alpha, 'g--')
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.2, alpha, 'r')
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.3, alpha, 'c')
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.4, alpha, 'b')
    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 0.5, alpha, 'm')

    read_spectrum_at_z('output/wab_test_7_t_0_nz_1281_nk_576.txt', 8, 2.0, alpha, 'k')

    k = np.logspace(-3,10,100)
    injection_slope = 4.66 - 1. / 3.
    k0 = 1. / 50.
    InverseGamma = 7e+5 * (k / k0)**(3.0 - injection_slope)
    #plt.plot(k, k**alpha * InverseGamma, 'k:')

    plt.xscale('log')
    plt.xlim([0.01,1e8])
    plt.xlabel(r'k [1/pc]',fontsize=28)

    plt.yscale('log')
    #plt.ylim([1e-8,1e8])
    plt.ylabel(r'$\Gamma^{-1}$(k) [Myr]',fontsize=28)

def plot_source_spectrum():
    alpha = 3.5
    filename = 'output/fcr_testloss_0_t_0_nz_3001_np_160.txt'
    read_spectrum_at_z(filename, 5, 0, alpha, 'r')
   
    plt.xscale('log')
    plt.xlim([1,1e7])
    plt.xlabel(r'p [GeV/c]',fontsize=30)
    
    plt.yscale('log')
    plt.ylabel(r'p$^3$ f(p) [m$^{-3}$ s$^{-1}$]',fontsize=30)

    return 'source_spectrum.pdf'

def plot_source_profile():
    filename = 'output/fcr_testloss_0_t_0_nz_3001_np_160.txt'
    read_profile_at_p(filename, 5, 1e2, 'r', False)
    read_profile_at_p(filename, 5, 1e3, 'r', False)
    read_profile_at_p(filename, 5, 1e4, 'r', False)

    plt.xscale('log')
    plt.xlabel(r'r [pc]',fontsize=30)

    plt.yscale('log')

    #plt.ylabel(r'f(p) [GeV$^{-3}$ m$^{-3}$ s$^{-1}$]',fontsize=30)

    return 'source_profile.pdf'

def plot_cascade_timescale():
    alpha = 0
    filename = 'output/fcr_geminga_5_t_3_nz_501_np_320.txt'
    
    read_spectrum_at_z(filename, 8,   5, alpha, 'y') # GROWTH
    read_spectrum_at_z(filename, 8,  10, alpha, 'g')
    read_spectrum_at_z(filename, 8,  20, alpha, 'r')

    read_spectrum_at_z(filename, 9,   5, alpha, 'y:') # DAMPING
    read_spectrum_at_z(filename, 9,  10, alpha, 'g:')
    read_spectrum_at_z(filename, 9,  20, alpha, 'r:')

    plt.plot([1, 1e6], [300, 300], 'k--')

    plt.legend(['5', '10', '20'], fontsize=20, frameon=False)

    plt.xscale('log')
    plt.xlim([10, 1e4])
    plt.xlabel(r'p [GeV/c]', fontsize=28)
    
    plt.yscale('log')
    plt.ylim([1e1, 1e5])
    plt.ylabel(r't [yr]', fontsize=28)

    return 'cascade_timescales.pdf'

def plot_timescales():
    alpha = 0
    
    plt.plot([0.1, 1e8], [300., 300.], 'k--') # age Geminga
    
    #    kpc = 3.1e21
    #vA = 10. * 1e5 # 10 km / s
    #t_ = (kpc / vA) / 3.1e10
    #plt.plot([0.1, 1e8], [t, t], 'r')
    
    filename = 'output/fcr_testloss_0_t_0_nz_3001_np_160.txt'
    read_spectrum_at_z(filename, 6, 0.0, alpha, 'g')
    read_spectrum_at_z(filename, 7, 0.0, alpha, 'b')
    
    plt.xscale('log')
    plt.xlim([1e2, 1e5])
    plt.xlabel(r'p [GeV/c]', fontsize=36)
    
    plt.yscale('log')
    plt.ylabel(r'timescale [kyr]', fontsize=36)

    return 'timescales.pdf'

def plot_positrons():
    alpha = 3
    plot_data('data/pos_FERMI.txt', alpha, 'c', 'FERMI-LAT', 1)
    plot_data('data/pos_PAMELA.txt', alpha, '#FFA500', 'PAMELA', 1)
    plot_data('data/pos_AMS02.txt', alpha, 'g', 'AMS-02', 1)

    filename = 'output/fcr_test_phase_2_t_340_nz_3001_np_160.txt'
    read_jcr_at_z(filename, 11, 250, 2, 'r')

    filename = 'output/fcr_test_phase_3_t_120_nz_3001_np_160.txt'
    read_jcr_at_z(filename, 11, 250, 2, 'b')
    
    plt.xscale('log')
    plt.xlim([1e1, 1e6])

    plt.yscale('log')
    plt.ylim([1e-2, 1e5])

    return 'positron_spectrum.pdf'

def plot_analytical_solution():
    alpha = 3.5

    filename = 'output/fcr_test_phase_3_t_0.2_nz_1001_np_160.txt'
    read_spectrum_at_z(filename, 4, .3, alpha, 'g')
    read_spectrum_at_z(filename, 5, .3, alpha, 'g:')

    filename = 'output/fcr_test_phase_3_t_0.4_nz_1001_np_160.txt'
    read_spectrum_at_z(filename, 4, .3, alpha, 'r')
    read_spectrum_at_z(filename, 5, .3, alpha, 'r:')
    
    filename = 'output/fcr_test_phase_3_t_0.6_nz_1001_np_160.txt'
    read_spectrum_at_z(filename, 4, .3, alpha, 'm')
    read_spectrum_at_z(filename, 5, .3, alpha, 'm:')

    #read_spectrum_at_z(filename, 10, 10,    alpha, 'g')
    #read_spectrum_at_z(filename, 11, 10, alpha, 'g:')

    plt.ylim([1e-10, 1e-1])

#filename = 'output/fcr_tested_3_t_5_nz_3001_np_160.txt'
#   read_spectrum_at_z(filename, 2, 10, alpha, 'r')

    plt.xscale('log')
    plt.yscale('log')

def plot_energy_density():
    filename = 'output/fcr_test_phase_9_t_324_nz_3001_np_160.txt'

    r = np.linspace(1, 100, 25)
    e = []
    for r_ in r:
        e_ = read_energy_density_at_z(filename, 2, r_)
        print r_, e_
        e.append(e_)

    plt.xlabel(r'z [pc]', fontsize=28)
    plt.ylabel(r'energy density [eV/cm$^3$]', fontsize=28)

    plt.xscale('log')
    plt.yscale('log')
    
    plt.plot(r, e)

print 'Plotting starts here...'

def plotit():
    
    #plt.figure(1)
    #plotname = plot_wave_spectrum()
    #plotname = plot_wave_profile()
    #plotname = plot_cascade_timescale()
    #plotname = plot_dzz_spectrum()
    #plotname = plot_dzz_profile()
    #plotname = plot_dzz_withtime()
    #plotname = plot_wave_profile()
    #plot_dzz_integral()
    #plot_dzz_timescale()
    #plot_giovanni()
    #plotname = plot_Jcr_spectrum()
    #plotname = plot_fcr_spectrum()
    plotname = plot_fcr_profile()
    #plotname = plot_dfdz_spectrum()
    #plotname = plot_dfdz_profile()
    #plotname = plot_flux_profile()
    #plot_Gamma_spectrum()
    #plotname = plot_source_spectrum()
    #plotname = plot_source_profile()
    #plotname = plot_timescales()
    #plotname = plot_positrons()
    #plotname = plot_analytical_solution()
    #plotname = plot_energy_density()

    #plt.show()
    plt.savefig(plotname, format='pdf', dpi=300)

plotit()
