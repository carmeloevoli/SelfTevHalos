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

    filename = 'output/fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128_t_0.txt'
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

def get_dzz_in_time(init_filename, t_max, energy):
    t = []
    D10, D20, D50 = [], [], []
    
    for i in range(1,t_max,1):
        t.append(i)
        filename = 'output/' + init_filename + '_t_' + str(i) + '.txt'
        #dzz = utils.read_spectrum_at_pz(filename, 3, 10, energy)
        z, dzz = np.loadtxt(filename, usecols=(0,3), skiprows=1, unpack=True)
        D10.append(dzz)
        #dzz = utils.read_spectrum_at_pz(filename, 3, 20, energy)
        D20.append(0.)
        #dzz = utils.read_spectrum_at_pz(filename, 3, 50, energy)
        D50.append(0.)
    
    txtfile = open(init_filename + '.txt', 'w+')
    size = len(t)
    for i in range(size):
        txtfile.write("%10.5e %10.5e %10.5e %10.5e\n" % (t[i], D10[i], D20[i], D50[i]))
    txtfile.close()
    #np.savez(init_filename, t=t, D10=D10, D20=D20, D50=D50)

def plot_dzz_withtime():
    fig = plt.figure(figsize=(10.5,8))
    ax = fig.add_subplot(111)
    
    get_dzz_in_time('fcr_test0.1pc_3.5_kol_3.8e33_nz_401_np_128', 600, 1e4)

    t, D10, D20, D50 = np.loadtxt('fcr_test0.1pc_3.5_kol_3.8e33_nz_401_np_128.txt', usecols=(0,1,2,3), unpack=True, skiprows=0)
    ax.plot(t, D10, color='tab:green', linestyle='--', label='source size = 0.1 pc')
    ax.plot(t, D20, color='tab:red', linestyle=':')
    ax.plot(t, D50, color='tab:green', linestyle=':')

    get_dzz_in_time('fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128', 200, 1e4)

    t, D10, D20, D50 = np.loadtxt('fcr_testN0.1pc_3.5_kol_3.8e33_nz_401_np_128.txt', usecols=(0,1,2,3), unpack=True, skiprows=0)
    ax.plot(t, D10, color='tab:blue', linestyle='--', label='source size = 1 pc')
    ax.plot(t, D20, color='tab:red', linestyle=':')
    ax.plot(t, D50, color='tab:green', linestyle=':')

    t, D = np.loadtxt('tim.txt', usecols=(0,1), unpack=True, skiprows=0)
    ax.plot(t, D, color='k', label='Tim')

    ax.set_xlabel(r't [kyr]')
    ax.set_xlim([1, 1e3])
    ax.set_xscale('log')

    ax.set_ylabel(r'D [cm$^2$/s]')
    ax.set_ylim([1e26, 1e30])
    ax.set_yscale('log')

    ax.legend()
    
    plt.savefig('dzz_with_time.pdf', format='pdf', dpi=300)

print ('Plotting starts here...')

if __name__ == "__main__":
#    plot_dzz_withtime()
    plot_dzz_spectrum()
#    plot_dzz_profile()
