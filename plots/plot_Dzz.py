import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_dzz_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    alpha = 0

    filename = 'output/fcr_test_3.5_kol_t_50_nz_601_np_128.txt'
    utils.read_spectrum_at_z(ax, filename, 3,   5, alpha, 'b')
    utils.read_spectrum_at_z(ax, filename, 3,  10, alpha, 'r')
    utils.read_spectrum_at_z(ax, filename, 3,  20, alpha, 'g')

    ax.legend(['5 pc', '10 pc', '20 pc'], fontsize=22, loc='upper left', frameon=False)

    ax.set_xlabel(r'p [GeV/c]')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])

    ax.set_ylabel(r'D(p) [cm$^2$/s]')
    ax.set_yscale('log')
    ax.set_ylim([1e26, 1e31])

    plt.savefig('dzz_spectrum.pdf', format='pdf', dpi=300)

def plot_dzz_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_test_3.5_kol_t_50_nz_601_np_128.txt'
    utils.read_profile_at_p(ax, filename, 3, 1e4, 'g', False)
    utils.read_profile_at_p(ax, filename, 3, 1e5, 'b', False)
    
    ax.legend([r'10 TeV', r'100 TeV'])

    ax.set_xlabel(r'z [pc]')
    ax.set_xscale('log')

    ax.set_ylabel(r'D(z) [cm$^2$/s]')
    ax.set_yscale('log')

    plt.savefig('dzz_profile.pdf', format='pdf', dpi=300)

def plot_dzz_withtime():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    t = []
    D10, D20, D50 = [], [], []
    
    for i in range(150):
        t.append(i)
        filename = 'output/fcr_test_3.5_kol_t_' + str(i) + '_nz_601_np_128.txt'
        dzz = utils.read_spectrum_at_pz(filename, 3, 10, 1e4)
        D10.append(dzz)
        dzz = utils.read_spectrum_at_pz(filename, 3, 20, 1e4)
        D20.append(dzz)
        dzz = utils.read_spectrum_at_pz(filename, 3, 50, 1e4)
        D50.append(dzz)
    
    ax.plot(t, D10, color='tab:blue')
    ax.plot(t, D20, color='tab:red')
    ax.plot(t, D50, color='tab:green')

    ax.set_xlabel(r't [kyr]')
    ax.set_xlim([0, 300])
    
    ax.set_ylabel(r'D(z) [cm$^2$/s]')
    ax.set_ylim([1e26, 1e30])
    ax.set_yscale('log')

    ax.legend(['10 pc', '20 pc', '50 pc'], fontsize=25, loc='lower right', frameon=False)

    plt.savefig('dzz_with_time.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__== "__main__":
    plot_dzz_spectrum()
    plot_dzz_profile()
    plot_dzz_withtime()
