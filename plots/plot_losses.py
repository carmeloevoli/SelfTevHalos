import matplotlib
matplotlib.use('MacOSX')
import matplotlib.pyplot as plt
plt.style.use('./my.mplstyle')
import utilities as utils
import numpy as np

def plot_losses():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
    utils.read_spectrum_at_z(ax, filename, 7, 0, 0, 'tab:red')

    ax.set_xlabel(r'p [GV]')
    ax.set_xscale('log')
    ax.set_xlim([1e2, 1e6])

    ax.set_ylabel(r'tau [kyr]')
    ax.set_yscale('log')
#    ax.set_ylim([1, 1e4])

    plt.savefig('losses.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__ == "__main__":
    plot_losses()
