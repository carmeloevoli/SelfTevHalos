import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_source_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    slope = 4
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
    utils.read_spectrum_at_z(ax, filename, 6,  0, slope, 'tab:red')

    ax.set_xlabel(r'p [GV]')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])

    ax.set_ylabel(r'p$^4$ Q(p) [GV cm$^{-3}$ s$^{-1}$]')
    ax.set_yscale('log')
    ax.set_ylim([1e-20, 1e-17])

    plt.savefig('source_spectrum.pdf', format='pdf', dpi=300)

def plot_source_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
    utils.read_profile_at_p(ax, filename, 6, 1e4, 'b', False)
    
    ax.set_xlabel(r'z [pc]')
    ax.set_xlim([-10, 10])
    
    ax.set_ylabel(r'Q(p) [a.u.]')

    plt.savefig('source_profile.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__ == "__main__":
    plot_source_spectrum()
    plot_source_profile()
