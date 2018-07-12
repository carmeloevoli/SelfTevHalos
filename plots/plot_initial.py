import matplotlib.pyplot as plt
import numpy as np
import plot_lib as pl

def plot_HAWC(pslope):
    x = 2e4
    y = 3.2e27 * (100. / 20.)**(1./3.)
    y_err_hi = 1.4e27 * (100. / 20.)**(1./3.)
    y_err_lo = 1.0e27 * (100. / 20.)**(1./3.)
    plt.errorbar(x, x**pslope * y, yerr=y_err_lo, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color='k')

def plot_source_spectrum():
    pslope = 4.0

    filename = 'output/fcr_pulsar_1D_t_0_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 7, 0, pslope, 'r')
    pl.read_spectrum_at_z(filename, 7, 0.1, pslope, 'b')
    pl.read_spectrum_at_z(filename, 7, 1, pslope, 'g')

    plt.xscale('log')
    plt.xlabel(r'p [GeV/c]',fontsize=30)
    
    plt.yscale('log')
    #plt.ylabel(r'p$^3$ f(p) []',fontsize=30)

def plot_source_profile():

    filename = 'output/fcr_pulsar_1D_t_0_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 7, 1e3, 'r', False)
    pl.read_profile_at_p(filename, 7, 1e4, 'b:', False)

    plt.xlim([0, 1])
    #plt.yscale('log')

def plot_dzz_spectrum():
    pslope = 0 # -1./3.

    filename = 'output/fcr_new_geminga_0_t_0_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'k:')
    
    filename = 'output/fcr_new_geminga_0_t_100_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'r')

    filename = 'output/fcr_new_geminga_3_t_100_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'b')

    filename = 'output/fcr_new_geminga_4_t_100_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'g')

    plt.legend([r'ISM', r'$\alpha$ = 3.5 - Kol', r'$\alpha$ = 3.2 - Kol', r'$\alpha$ = 3.5 - Kra'], loc='upper left', fontsize=18)
    
    filename = 'output/fcr_new_geminga_0_t_50_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'r--')
    
    filename = 'output/fcr_new_geminga_3_t_50_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'b--')
    
    filename = 'output/fcr_new_geminga_4_t_50_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'g--')

    plt.text(1.5e3, 1.05e27, 't = 50 kyr', fontsize=21, rotation=20)
    plt.text(.2e3, 2.5e27, 't = 100 kyr', fontsize=21, rotation=25)
    
    plt.xscale('log')
    plt.xlabel(r'p [GeV/c]',fontsize=30)
    
    plt.yscale('log')
    plt.ylabel(r'D [cm$^2$/s]',fontsize=30)
    plt.ylim([1e26, 1e31])
    
    return 'dzz_spectrum.pdf'

def plot_dzz_profile():
    filename = 'output/fcr_new_geminga_0_t_100_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e2, 'r', False)
    pl.read_profile_at_p(filename, 3, 1e4, 'r:', False)

    filename = 'output/fcr_new_geminga_3_t_100_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e2, 'b', False)
    pl.read_profile_at_p(filename, 3, 1e4, 'b:', False)

    filename = 'output/fcr_new_geminga_4_t_100_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e2, 'g', False)
    pl.read_profile_at_p(filename, 3, 1e4, 'g:', False)

    plt.xlabel(r'z [pc]',fontsize=30)

    plt.ylabel(r'D [cm2/s]',fontsize=30)
    plt.yscale('log')

    return 'dzz_profile_2.pdf'
    
def plot_dzz_time():
    pslope = 0

    plot_HAWC(pslope)

    filename = 'output/fcr_geminga_5_t_0_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'k:')
    filename = 'output/fcr_geminga_5_t_2_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'y')
    filename = 'output/fcr_geminga_5_t_4_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'g')
    filename = 'output/fcr_geminga_5_t_6_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'r')
    filename = 'output/fcr_geminga_5_t_8_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'b')
    filename = 'output/fcr_geminga_5_t_10_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 3, 10, pslope, 'c')   
    #filename = 'output/fcr_geminga_5_t_10_nz_1201_np_256.txt'
    #pl.read_spectrum_at_z(filename, 3, 10, pslope, 'm')   

    plt.legend(['ISM','10','20','30','40','50', '90'])

    plt.xscale('log')
    plt.xlabel(r'p [GeV/c]')
    
    plt.yscale('log')
    plt.ylabel(r'D [cm$^2$/s]')

    return 'dzz_time_Kol_vs_Kra.pdf'

def plot_dzz_profile_time():
    filename = 'output/fcr_geminga_5_t_0_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'k:', False)
    filename = 'output/fcr_geminga_5_t_2_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'y', False)
    filename = 'output/fcr_geminga_5_t_4_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'g', False)
    filename = 'output/fcr_geminga_5_t_6_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'r', False)
    filename = 'output/fcr_geminga_5_t_8_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'b', False)
    filename = 'output/fcr_geminga_5_t_10_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 3, 1e4, 'c', False)

#    pl.read_profile_at_p(filename, 3, 1e4, 'm', False)
    
    plt.legend(['ISM','5','10','20','30','40','50'])
    
    plt.xlabel(r'z [pc]')
    plt.xlim([0,5])
    
    plt.yscale('log')
    plt.ylabel(r'D [cm$^2$/s]')

    return 'dzz_profile_time.pdf'

def plot_fcr_time():
    pslope = 3

    filename = 'output/fcr_new_geminga_0_t_100_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 2, 10, pslope, 'y')
    filename = 'output/fcr_new_geminga_0_t_200_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 2, 10, pslope, 'g')
    filename = 'output/fcr_new_geminga_0_t_300_nz_1201_np_256.txt'
    pl.read_spectrum_at_z(filename, 2, 10, pslope, 'r')
    
    plt.legend(['100','200','300'])

    #pl.plot_data('AMS02.txt', 1, 'g', 'AMS-02')
    
    plt.xscale('log')
    plt.xlabel(r'p [GeV/c]')
    plt.xlim([1e2, 1e5])
    
    plt.yscale('log')
    plt.ylabel(r'p$^3$ f []')
    plt.ylim([1e-11, 1e-5])

    return 'fcr_time.pdf'

def plot_fcr_profile_time():
    filename = 'output/fcr_pulsar_1D_t_0_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'k:', False)
    filename = 'output/fcr_pulsar_1D_t_5_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'y', False)
    filename = 'output/fcr_pulsar_1D_t_10_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'g', False)
    filename = 'output/fcr_pulsar_1D_t_20_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'r', False)
    filename = 'output/fcr_pulsar_1D_t_30_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'b', False)
    filename = 'output/fcr_pulsar_1D_t_40_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'c', False)
    filename = 'output/fcr_pulsar_1D_t_50_nz_1201_np_256.txt'
    pl.read_profile_at_p(filename, 2, 1e4, 'm', False)
    
    plt.legend(['ISM','5','10','20','30','40','50'])
    
    plt.xlabel(r'z [pc]')
    
    plt.yscale('log')
    plt.ylabel(r'f []')
    plt.ylim([1e-21, 1e-17])
    
    return 'fcr_profile_time.pdf'

def calc_dzz_in_time(energy, distance, filename_id):
    t = []
    f = open(filename_id + '_' + str(energy) + '_' + str(distance) + '_dzz_time.txt', 'wb')
    for i in range(300):
        print i
        t.append(i)
        filename = 'output/fcr_' + filename_id + '_t_' + str(i) + '_nz_1201_np_256.txt'
        p, D = pl.get_spectrum_at_z(filename, 3, distance)
        logDzz_1 = np.interp(np.log(1e2), np.log(p), np.log(D))
        logDzz_2 = np.interp(np.log(1e3), np.log(p), np.log(D))
        logDzz_3 = np.interp(np.log(1e4), np.log(p), np.log(D))
        logDzz_4 = np.interp(np.log(1e5), np.log(p), np.log(D))
        f.write("%i  %.4e  %.4e  %.4e  %.4e\n" % (i, np.exp(logDzz_1), np.exp(logDzz_2), np.exp(logDzz_3), np.exp(logDzz_4)))
    f.close()

def plot_dzz_in_time_new():
    id_energy = 3

    t, Dzz = np.loadtxt('new_geminga_5_0_20_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)
    plt.plot(t, Dzz, color='b', label='10 pc')
    
    t, Dzz_min = np.loadtxt('new_geminga_6_0_20_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)
    t, Dzz_max = np.loadtxt('new_geminga_7_0_20_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)

    plt.fill_between(t, Dzz_min, Dzz_max, facecolor='dodgerblue', alpha=0.5, linewidth=0)
    
    t, Dzz = np.loadtxt('new_geminga_5_0_10_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)
    plt.plot(t, Dzz, color='r', label='20 pc')

    plt.plot(t, np.max(Dzz) * Dzz / Dzz, 'k:', label='ISM')

    t, Dzz_min = np.loadtxt('new_geminga_6_0_10_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)
    t, Dzz_max = np.loadtxt('new_geminga_7_0_10_dzz_time.txt', skiprows=0, usecols=(0,id_energy), unpack=True)

    plt.fill_between(t, Dzz_min, Dzz_max, facecolor='tomato', alpha=0.5, linewidth=0)
    
    plt.text(200, 2.0e26, r'$\alpha$ = 3.2', fontsize=25)
    #plt.text(200, 2.0e26, r'$\alpha$ = 3.5', fontsize=25)

    plt.text(200, 3.5e26, 'Kraichnan', fontsize=25)
    #plt.text(200, 3.5e26, 'Kolmogorov', fontsize=25)

    #plt.text(200, 4.0e29, 'E = 100 GeV', fontsize=25)
    plt.text(200, 2.0e29, 'E = 10 TeV', fontsize=25)

    #plt.text(50,  2e29, 'ISM', fontsize=25)
    #plt.text(50,  0.45e30, 'ISM', fontsize=25)

    #plt.text(250, 1e28, r'10 pc', color='r', fontsize=22)
    #plt.text(230, .4e28, r'20 pc', color='b', fontsize=22)
    #plt.text(250, 1.9e29, r'10 pc', color='r', fontsize=22)
    #plt.text(230, .45e29, r'20 pc', color='b', fontsize=22)

    plt.legend(loc='lower left')
    
    plt.yscale('log')
    plt.ylabel(r'D [cm$^2$/s]')
    plt.ylim([1e26, 1e30])

    plt.xlabel(r't [kyr]')
    
    return 'dzz_time_3_2_kra_1e4.pdf'
    #return 'dzz_time_3_5_kol_1e2.pdf'

def plot_dzz_in_time(energy):
    t = []
    
    Dzz_0 = []
    Dzz_1 = []
    Dzz_2 = []

    Dzz_3 = []
    Dzz_4 = []
    Dzz_5 = []
       
    for i in range(300):
        t.append(i)
        
        filename = 'output/fcr_new_geminga_6_t_' + str(i) + '_nz_1201_np_256.txt'
        p, D = pl.get_spectrum_at_z(filename, 3, 10)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_0.append(np.exp(logDzz))
        p, D = pl.get_spectrum_at_z(filename, 3, 20)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_3.append(np.exp(logDzz))
        
        filename = 'output/fcr_new_geminga_5_t_' + str(i) + '_nz_1201_np_256.txt'
        p, D = pl.get_spectrum_at_z(filename, 3, 10)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_1.append(np.exp(logDzz))
        p, D = pl.get_spectrum_at_z(filename, 3, 20)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_4.append(np.exp(logDzz))
        
        filename = 'output/fcr_new_geminga_7_t_' + str(i) + '_nz_1201_np_256.txt'
        p, D = pl.get_spectrum_at_z(filename, 3, 10)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_2.append(np.exp(logDzz))
        p, D = pl.get_spectrum_at_z(filename, 3, 20)
        logDzz = np.interp(np.log(energy), np.log(p), np.log(D))
        Dzz_5.append(np.exp(logDzz))
        
        print (i)
        
    #plt.plot(t, Dzz_0, color='r', label=r'B$_0$ = 0.5 $\mu$G')
    #plt.plot(t, Dzz_1, color='b', label=r'B$_0$ = 1 $\mu$G')
    #plt.plot(t, Dzz_2, color='g', label=r'B$_0$ = 2 $\mu$G')

    plt.plot(t, Dzz_0, color='r', label=r'$\alpha$ = 3.2')
    plt.plot(t, Dzz_1, color='b', label=r'fiducial')
    plt.plot(t, Dzz_2, color='g', label=r'Kraichnan')
    
    plt.plot(t, Dzz_3, color='r', linestyle=':')
    plt.plot(t, Dzz_4, color='b', linestyle=':')
    plt.plot(t, Dzz_5, color='g', linestyle=':')

    #plt.legend()

    #plt.title(r'$\alpha$ = 3.2 - 
    
    plt.yscale('log')
    plt.ylabel(r'D [cm$^2$/s]')
    plt.ylim([1e26, 1e30])

    plt.xlabel(r't [kyr]')
    
    return 'dzz_time_3.2_kra_1e2.pdf'

pl.set_plot_style()
    
#plot_source_spectrum()
#plot_source_profile()
filename = plot_dzz_spectrum()
#filename = plot_dzz_profile()
#filename = plot_dzz_time()
#filename = plot_dzz_profile_time()
#filename = plot_fcr_time()
#filename = plot_fcr_profile_time()
#filename = plot_dzz_in_time_new()

#calc_dzz_in_time(0, 10, 'new_geminga_0')
#calc_dzz_in_time(0, 20, 'new_geminga_0')
#calc_dzz_in_time(0, 10, 'new_geminga_1')
#calc_dzz_in_time(0, 20, 'new_geminga_1')
#calc_dzz_in_time(0, 10, 'new_geminga_2')
#calc_dzz_in_time(0, 20, 'new_geminga_2')

#plt.show()
plt.savefig(filename)
