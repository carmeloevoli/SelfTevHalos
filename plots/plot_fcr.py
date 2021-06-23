import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_fcr_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    alpha = 3
    p = np.logspace(-2, 6, 100)
    
    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_spectrum_at_z(ax, filename, 2,  5, alpha, 'y')
    utils.read_spectrum_at_z(ax, filename, 2, 10, alpha, 'g')
    utils.read_spectrum_at_z(ax, filename, 2, 50, alpha, 'r')

    ax.legend(['z=5', 'z=10', 'z=50'])

#    filename = 'output/fcr_test_3D_0_t_30_nz_401_np_20.txt'
#    read_spectrum_at_z(filename, 2,    5, alpha, 'y:')
#    read_spectrum_at_z(filename, 2,   10, alpha, 'g:')
#    read_spectrum_at_z(filename, 2,   50, alpha, 'r:')

    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])
    ax.set_xlabel(r'p [GeV/c]',fontsize=30)
    
    ax.set_yscale('log')
    ax.set_ylim([1e-10, 1e-4])
    ax.set_ylabel(r'p$^3$ f(p) [m$^{-3}$]',fontsize=30)

    p = np.logspace(2, 5, 100)
    #plt.plot(p, 0.5e-11 * (p / 1e4)**(-3.5), 'k:')
    #plt.plot(p, 1.0e-7 * (p / 1e2)**(-1.5), 'k:')
    plt.savefig('fcr_spectrum.pdf', format='pdf', dpi=300)

def plot_fcr_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_profile_at_p(ax, filename, 2, 1e5, 'r', True)
    utils.read_profile_at_p(ax, filename, 2, 1e2, 'm', True)

#    filename = 'output/fcr_test_3D_0_t_30_nz_401_np_20.txt'
#    read_profile_at_p(filename, 2, 1e5, 'r:', True)
#    read_profile_at_p(filename, 2, 1e2, 'm:', True)
    
#    filename = 'output/fcr_test_1D_1_t_4_nz_1001_np_20.txt'
#    read_profile_at_p(filename, 2, 1e5, 'r:', True)
#    read_profile_at_p(filename, 2, 1e3, 'b:', True)
    
    #plt.xscale('log')
    #plt.xlim([0, 10])
    ax.set_xlabel(r'z [pc]', fontsize=28)
    
    #plt.yscale('log')
    #plt.ylim([0.8, 1])
    ax.set_ylabel(r'f(z) [au]', fontsize=28)
    
    #z = np.logspace(-1, 3, 100)
    #plt.plot(z, np.exp(-0.5 * (z / 25.)**2), 'k--')

#plt.legend(['10 GeV', '100 GeV', '1 TeV', '10 TeV'], fontsize=25, loc='upper right', frameon=False)

    plt.savefig('fcr_profile.pdf', format='pdf', dpi=300)

def plot_dzz_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    alpha = 0 # -1. / 3.

    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_spectrum_at_z(ax, filename, 3,   5, alpha, 'b')
    utils.read_spectrum_at_z(ax, filename, 3,  10, alpha, 'r')
    utils.read_spectrum_at_z(ax, filename, 3,  20, alpha, 'g')

#    filename = 'output/fcr_test_1D_0_t_300_nz_1201_np_96.txt'
#    read_spectrum_at_z(filename, 3,  5, alpha, 'b:')
#    read_spectrum_at_z(filename, 3, 10, alpha, 'r:')
#    read_spectrum_at_z(filename, 3, 20, alpha, 'g:')

    ax.legend(['5 pc', '10 pc', '20 pc'], fontsize=22, loc='upper left', frameon=False)

#    #filename = 'output/fcr_test_1D_0_t_0_nz_1201_np_96.txt'
#    #read_spectrum_at_z(filename, 3, 10, alpha, 'k:')
#
#    x = 2e4
#    y = 3.2e27 * (100. / 20.)**(1./3.)
#    y_err_hi = 1.4e27 * (100. / 20.)**(1./3.)
#    y_err_lo = 1.0e27 * (100. / 20.)**(1./3.)
#
#    plt.errorbar(x, x**alpha * y, yerr=y_err_lo, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color='k')

    #plt.plot([5e4, 5e4], [1e20, 1e40], 'c:')

#    p = np.logspace(3,6,100)
#
#    c = 3e10 # cm / s
#    rL = 3.1e18 * (p / 1e6)
#
#    plt.plot(p, rL * c / 3.0, 'k--')
    
    #plt.text(1e4, 2e27, 'H', fontsize=20)
    
    #ax.text(1e3, 2e30, 'Galactic Diffusion', fontsize=25, rotation=11)
    #ax.text(1e4, 5e26, 'Bohm limit', fontsize=25, rotation=30)

    ax.set_xscale('log')
    #plt.xlim([1e1, 1e6])
    ax.set_xlabel(r'p [GeV/c]',fontsize=35,labelpad=10)

    ax.set_yscale('log')
    #plt.ylim([1e26, 1e31])
    ax.set_ylabel(r'D(p) [cm$^2$/s]',fontsize=35,labelpad=15)

    plt.savefig('dzz_spectrum.pdf', format='pdf', dpi=300)

def plot_dzz_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_profile_at_p(ax, filename, 3, 1e5, 'b', False)
    utils.read_profile_at_p(ax, filename, 3, 1e3, 'g', False)
    
    ax.legend([r'$1e5$', r'$1e3$'])
    #read_profile_at_p(filename, 3, 1e4, 'r', True)
    #read_profile_at_p(filename, 3, 1e5, 'b', True)

#    filename = 'output/fcr_test_1D_1_t_20_nz_1201_np_96.txt'
#    read_profile_at_p(filename, 3, 1e4, 'b:', False)

    ax.set_xlabel(r'z [pc]', fontsize=38)
    #plt.xlim([0, 10])
    ax.set_xscale('log')

    ax.set_ylabel(r'D(z) [cm$^2$/s]', fontsize=38)
    #plt.ylim([0.01,1.1])
    ax.set_yscale('log')

#    ax.legend([r'$E = 100$ TeV'], fontsize=25, loc='upper right', frameon=False)
    plt.savefig('dzz_profile.pdf', format='pdf', dpi=300)

def plot_dfdz_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'

    alpha = 3.0
    utils.read_spectrum_at_z(ax, filename, 4,   5, alpha, 'b')
    utils.read_spectrum_at_z(ax, filename, 4,  10, alpha, 'r')
    utils.read_spectrum_at_z(ax, filename, 4,  50, alpha, 'g')
    utils.read_spectrum_at_z(ax, filename, 4, 100, alpha, 'm')

    ax.legend(['5 pc', '10 pc', '50 pc', '100 pc'], fontsize=22, loc='best', frameon=False)

#    read_spectrum_at_z(filename, 4, -1, alpha, 'g')
#    read_spectrum_at_z(filename, 4, -2, alpha, 'c')
#    read_spectrum_at_z(filename, 4, -5, alpha, 'b')
#    read_spectrum_at_z(filename, 4,  0, alpha, 'y')
#    read_spectrum_at_z(filename, 4,  1, alpha, 'g:')
#    read_spectrum_at_z(filename, 4,  2, alpha, 'c:')
#    read_spectrum_at_z(filename, 4,  5, alpha, 'b:')

    plt.xscale('log')
    #plt.xlim([0.1,1e7])
    plt.xlabel(r'p [GeV/c]',fontsize=28)
    
    plt.yscale('log')
    plt.ylim([1e-10,1e-2])
    plt.ylabel(r'p3 df/dz [m$^{-3}$ kpc$^{-1}$]',fontsize=28)
    plt.savefig('dfdz_spectrum.pdf', format='pdf', dpi=300)

def plot_dfdz_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'

    utils.read_profile_at_p(ax, filename, 4, 1e4, 'r', False)
    utils.read_profile_at_p(ax, filename, 4, 1e5, 'b', False)

    #read_profile_at_p(filename, 4, 1e2, 'g', True)
    #read_profile_at_p(filename, 4, 1e3, 'c', True)
    #read_profile_at_p(filename, 4, 1e4, 'b', True)
    #read_profile_at_p(filename, 4, 1e5, 'm', False)

#    filename = 'output/fcr_test_3D_1_t_1_nz_401_np_20.txt'
#    read_profile_at_p(filename, 4, 1e5, 'b', False)
#
#    filename = 'output/fcr_test_3D_1_t_2_nz_401_np_20.txt'
#    read_profile_at_p(filename, 4, 1e5, 'r', False)
    
    ax.legend([r'$10^4$ GeV', '$10^5$ GeV'], fontsize=22, loc='best', frameon=False)

    #plt.xscale('log')
    plt.xlim([0, 25])
    plt.xlabel(r'z [kpc]',fontsize=28)
    
    plt.yscale('log')
    #plt.ylim([1e-5,1e1])
    plt.ylabel(r'p$^3$ df/dz [m$^{-3}$ kpc$^{-1}$]',fontsize=28)

    plt.savefig('dfdz_profile.pdf', format='pdf', dpi=300)

def plot_waves_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'

    alpha = 0.0
    utils.read_spectrum_at_z(ax, filename, 6,   5, alpha, 'b')
    utils.read_spectrum_at_z(ax, filename, 6,  10, alpha, 'r')
    utils.read_spectrum_at_z(ax, filename, 6,  50, alpha, 'g')
    utils.read_spectrum_at_z(ax, filename, 6, 100, alpha, 'm')

    ax.legend(['5 pc', '10 pc', '50 pc', '100 pc'], fontsize=22, loc='best', frameon=False)

#    read_spectrum_at_z(filename, 4, -1, alpha, 'g')
#    read_spectrum_at_z(filename, 4, -2, alpha, 'c')
#    read_spectrum_at_z(filename, 4, -5, alpha, 'b')
#    read_spectrum_at_z(filename, 4,  0, alpha, 'y')
#    read_spectrum_at_z(filename, 4,  1, alpha, 'g:')
#    read_spectrum_at_z(filename, 4,  2, alpha, 'c:')
#    read_spectrum_at_z(filename, 4,  5, alpha, 'b:')

    plt.xscale('log')
    #plt.xlim([0.1,1e7])
    plt.xlabel(r'p [GeV/c]',fontsize=28)
    
    plt.yscale('log')
    #plt.ylim([1e-10,1e-2])
    plt.ylabel(r'W [m]',fontsize=28)
    plt.savefig('waves_spectrum.pdf', format='pdf', dpi=300)

def plot_waves_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'

    utils.read_profile_at_p(ax, filename, 6, 1e4, 'r', False)
    utils.read_profile_at_p(ax, filename, 6, 1e5, 'b', False)

    #read_profile_at_p(filename, 4, 1e2, 'g', True)
    #read_profile_at_p(filename, 4, 1e3, 'c', True)
    #read_profile_at_p(filename, 4, 1e4, 'b', True)
    #read_profile_at_p(filename, 4, 1e5, 'm', False)

#    filename = 'output/fcr_test_3D_1_t_1_nz_401_np_20.txt'
#    read_profile_at_p(filename, 4, 1e5, 'b', False)
#
#    filename = 'output/fcr_test_3D_1_t_2_nz_401_np_20.txt'
#    read_profile_at_p(filename, 4, 1e5, 'r', False)
    
    ax.legend([r'$10^4$ GeV', '$10^5$ GeV'], fontsize=22, loc='best', frameon=False)

    #plt.xscale('log')
    plt.xlim([0, 25])
    plt.xlabel(r'z [kpc]',fontsize=28)
    
    plt.yscale('log')
    #plt.ylim([1e-5,1e1])
    plt.ylabel(r'W [m]',fontsize=28)

    plt.savefig('waves_profile.pdf', format='pdf', dpi=300)



    





    
    
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





def plot_flux_profile():
    filename = 'output/fcr_test_1D_0_t_10_nz_1001_np_20.txt'
    read_profile_at_p(filename, 5, 1e2, 'g', False)
    read_profile_at_p(filename, 5, 1e3, 'c', False)
    read_profile_at_p(filename, 5, 1e4, 'b', False)
    read_profile_at_p(filename, 5, 1e5, 'm', False)
    
    plt.xlim([90, 100])
    plt.yscale('log')


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

def plot_electrons():
    alpha = 3
    plot_data('data/ele_AMS02.txt', alpha, 'g', 'AMS-02', 1)
    
    f0 = 2e-1
    p0 = 10.
    gamma = 3.2
    p = np.logspace(1, 5)
    f = f0 * (p / p0)**(-gamma)

    plt.plot(p, p**alpha * f)

    plt.xscale('log')
    plt.xlim([1e1, 1e5])
    
    plt.yscale('log')
    plt.ylim([1e1, 1e3])
    
    return 'electron_spectrum.pdf'

def plot_analytical_solution():

    filename = 'output/test_test_3D_2_t_0.1_nz_401_np_20.txt'
    read_profile_at_p(filename, 2, 1e2, 'g', False)
    read_profile_at_p(filename, 3, 1e2, 'r:', False)

#read_spectrum_at_z(filename, 2, 10, alpha, 'r')
#    read_spectrum_at_z(filename, 3, 10, alpha, 'r:')
    
    #    read_spectrum_at_z(filename, 2, 30, alpha, 'm')
    #read_spectrum_at_z(filename, 3, 30, alpha, 'm:')

    #read_spectrum_at_z(filename, 10, 10,    alpha, 'g')
    #read_spectrum_at_z(filename, 11, 10, alpha, 'g:')

    #plt.ylim([1e-10, 1e-1])

#filename = 'output/fcr_tested_3_t_5_nz_3001_np_160.txt'
#   read_spectrum_at_z(filename, 2, 10, alpha, 'r')

    plt.xscale('log')
    plt.yscale('log')

    return 'analytical_solution.pdf'

def plot_energy_density():
    filename = 'output/fcr_test_phase_9_t_324_nz_3001_np_160.txt'

    r = np.linspace(1, 100, 25)
    e = []
    for r_ in r:
        e_ = read_energy_density_at_z(filename, 2, r_)
        print (r_, e_)
        e.append(e_)

    plt.xlabel(r'z [pc]', fontsize=28)
    plt.ylabel(r'energy density [eV/cm$^3$]', fontsize=28)

    plt.xscale('log')
    plt.yscale('log')
    
    plt.plot(r, e)

print ('Plotting starts here...')

if __name__== "__main__":
    plot_fcr_spectrum()
    plot_fcr_profile()
    plot_dzz_spectrum()
    plot_dzz_profile()
    plot_dfdz_spectrum()
    plot_dfdz_profile()
    plot_waves_spectrum()
    plot_waves_profile()

    #plot_energy_density()
    #plotname = plot_cascade_timescale()
    #plotname = plot_dzz_withtime()
    #plotname = plot_wave_profile()
    #plot_dzz_integral()
    #plot_dzz_timescale()
    #plot_giovanni()
    #plotname = plot_Jcr_spectrum()

    #plotname = plot_flux_profile()
    #plot_Gamma_spectrum()
    #plotname = plot_source_spectrum()
    #plotname = plot_source_profile()
    #plotname = plot_timescales()
    #plotname = plot_positrons()
    #plotname = plot_electrons()
    #plotname = plot_analytical_solution()

