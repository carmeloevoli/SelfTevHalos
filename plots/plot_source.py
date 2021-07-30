import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_source_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    alpha = 4

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
    utils.read_spectrum_at_z(ax, filename, 7,  0, alpha, 'tab:red')

    #ax.legend(['0 pc', '1 pc'], fontsize=22, loc='upper left', frameon=False)

    ax.set_xlabel(r'p [GeV/c]')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])

    ax.set_ylabel(r'p$^4$ Q(p) [a.u.]')
    ax.set_yscale('log')
    ax.set_ylim([1e42, 1e44])

    plt.savefig('source_spectrum.pdf', format='pdf', dpi=300)

def plot_source_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
    utils.read_profile_at_p(ax, filename, 7, 1e4, 'b', False)
    
    #ax.legend([r'10 TeV', r'100 TeV'])

    ax.set_xlabel(r'z [pc]')
    ax.set_ylabel(r'Q(p) [a.u.]')

    plt.savefig('source_profile.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__ == "__main__":
    plot_source_spectrum()
    plot_source_profile()
