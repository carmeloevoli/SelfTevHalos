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

    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_spectrum_at_z(ax, filename, 3,   5, alpha, 'b')
    utils.read_spectrum_at_z(ax, filename, 3,  10, alpha, 'r')
    utils.read_spectrum_at_z(ax, filename, 3,  20, alpha, 'g')

    ax.legend(['5 pc', '10 pc', '20 pc'], fontsize=22, loc='upper left', frameon=False)

    ax.set_xscale('log')
    ax.set_xlabel(r'p [GeV/c]',fontsize=35,labelpad=10)

    ax.set_yscale('log')
    ax.set_ylabel(r'D(p) [cm$^2$/s]',fontsize=35,labelpad=15)

    plt.savefig('dzz_spectrum.pdf', format='pdf', dpi=300)

def plot_dzz_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_test_3.5_kol_t_5_nz_1201_np_256.txt'
    utils.read_profile_at_p(ax, filename, 3, 1e3, 'g', False)
    utils.read_profile_at_p(ax, filename, 3, 1e5, 'b', False)
    
    ax.legend([r'1 TeV', r'100 TeV'])

    ax.set_xlabel(r'z [pc]', fontsize=38)
    ax.set_xscale('log')

    ax.set_ylabel(r'D(z) [cm$^2$/s]', fontsize=38)
    ax.set_yscale('log')

    plt.savefig('dzz_profile.pdf', format='pdf', dpi=300)

def plot_dzz_withtime():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    t = []
    D0, D1, D2 = [], [], []
    
    for i in range(20):
        t.append(i)
        filename = 'output/fcr_test_3.5_kol_t_' + str(i) + '_nz_1201_np_256.txt'
        dzz = utils.read_spectrum_at_pz(filename, 3, 0.01, 1e4)
        D0.append(dzz)
        dzz = utils.read_spectrum_at_pz(filename, 3, 0.05, 1e4)
        D1.append(dzz)
        dzz = utils.read_spectrum_at_pz(filename, 3, 0.25, 1e4)
        D2.append(dzz)
    
    ax.plot(t, D0, color='tab:green')
    ax.plot(t, D1, color='tab:red')
    ax.plot(t, D2, color='tab:blue')

    ax.set_xlabel(r't [kyr]')
    ax.set_xlim([0,100])
    
    ax.set_ylabel(r'D(z) [cm$^2$/s]')
    ax.set_ylim([1e26,1e30])
    ax.set_yscale('log')

    ax.legend(['10 pc', '50 pc', '250 pc'], fontsize=25, loc='upper right', frameon=False)

    plt.savefig('dzz_with_time.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__== "__main__":
#    plot_dzz_spectrum()
#    plot_dzz_profile()
    plot_dzz_withtime()
