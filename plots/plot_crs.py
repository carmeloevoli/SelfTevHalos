import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_fcr_spectrum():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    slope = 3.5
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_1.txt'
    utils.read_spectrum_at_z(ax, filename, 2,  20, slope, 'tab:blue')

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_10.txt'
    utils.read_spectrum_at_z(ax, filename, 2,  20, slope, 'tab:orange')
    
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_50.txt'
    utils.read_spectrum_at_z(ax, filename, 2,  20, slope, 'tab:red')

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_100.txt'
    utils.read_spectrum_at_z(ax, filename, 2,  20, slope, 'tab:purple')
    
    p = np.logspace(2, 6, 100)
    ax.plot(p, 1.2e-10 * np.power(p / 1e2, -1./3.), ':')
    
    ax.set_xlabel(r'p [GV]')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])

    ax.set_ylabel(r'p$^4$ f(p) [GV cm$^{-3}$]')
    ax.set_yscale('log')
    ax.set_ylim([1e-12, 3e-10])

    plt.savefig('fcr_spectrum.pdf', format='pdf', dpi=300)

def plot_fcr_profile():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_10.txt'
    utils.read_profile_at_p(ax, filename, 2, 1e3, 'tab:blue', True)

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_50.txt'
    utils.read_profile_at_p(ax, filename, 2, 1e3, 'tab:orange', True)

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_100.txt'
    utils.read_profile_at_p(ax, filename, 2, 1e3, 'tab:red', True)
    
    ax.set_xlabel(r'z [pc]')
    ax.set_xlim([-100, 100])
    
    ax.set_ylabel(r'f(p) [a.u.]')
    ax.set_ylim([0, 1.1])
    
    plt.savefig('fcr_profile.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__ == "__main__":
    plot_fcr_spectrum()
    plot_fcr_profile()
