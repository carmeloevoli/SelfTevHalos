import numpy as np

def plot_data(filename, slope_plot, color, label, norm):
    T, y, err_y_lo, err_y_up = np.loadtxt(filename,skiprows=2,usecols=(0,3,8,9),unpack=True)
    y_err = [T**slope_plot * err_y_lo, T**slope_plot * err_y_up]
    plt.errorbar(T, T**slope_plot * y, yerr=y_err, fmt='o', markersize='8', elinewidth=2, capsize=6, capthick=2, color=color)


def read_spectrum_at_pz(filename, column, z_search, p_search):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (p_search < min(p) or p_search > max(p)):
        print ("WARNING: p is out of range!")
    if (z_search < min(z) or z_search > max(z)):
        print ("WARNING: p is out of range!")
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
    print ('f is',f)
    return f
  
  
def read_profile_at_p(ax, filename, column, p_search, linestyle, donormalize):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (p_search < min(p) or p_search > max(p)):
        print ("WARNING: p is out of range!")
    best_p = -1000.0
    for i in range(len(p)):
        if abs(p[i] - p_search) < abs(best_p - p_search):
            best_p = p[i]
    ind = (p == best_p).nonzero()
    z = z[ind]
    f = f[ind]
    print ("plot f in the range = ", min(f), max(f), np.mean(f))
    if donormalize:
        ax.plot(z, f / max(f), linestyle)
    else:
        ax.plot(z, f, linestyle)


def read_spectrum_at_z(ax, filename, column, z_search, alpha, linestyle):
    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
    if (z_search < min(z) or z_search > max(z)):
        print ("WARNING: z is out of range!")
    best_z = -1000.0
    for i in range(len(z)):
        if abs(z[i] - z_search) < abs(best_z - z_search):
            best_z = z[i]
    ind = (z == best_z).nonzero()
    p = p[ind] # GeV/c
    f = f[ind] # GeV^-1 m^-2 s^-1
    print ("plot f in the range = ", min(f), max(f), np.mean(f))
    ax.plot(p, p**alpha * f, linestyle)
    
    
#def read_energy_density_at_z(filename, column, z_search):
#    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
#    if (z_search < min(z) or z_search > max(z)):
#        print ("WARNING: z is out of range!")
#    best_z = -1000.0
#    for i in range(len(z)):
#        if abs(z[i] - z_search) < abs(best_z - z_search):
#            best_z = z[i]
#    ind = (z == best_z).nonzero()
#    p = p[ind] # GeV/c
#    f = f[ind] # GeV/c^-3 m^-3
#
#    epsilon = np.sum(p**4 * f)
#    dlnp = np.log(p[1] / p[0])
#
#    return 4.0 * 3.14 * epsilon * dlnp * 1e3


#def read_v_at_z(filename, z_search, alpha, linestyle):
#    z, p, f, D, dfdz = np.loadtxt(filename, skiprows=1, usecols=(0,1,2,3,4), unpack=True)
#    if (z_search < min(z) or z_search > max(z)):
#        print ("WARNING: z is out of range!")
#    best_z = -1000.0
#    for i in range(len(z)):
#        if abs(z[i] - z_search) < abs(best_z - z_search):
#            best_z = z[i]
#    ind = (z == best_z).nonzero()
#    p = p[ind] # GeV/c
#    f = f[ind] # GeV^-1 m^-2 s^-1
#    D = D[ind]
#    dfdz = dfdz[ind]
#    plt.plot(p, p**alpha * D / f * dfdz / 3e21, linestyle)


#def read_timescale(filename, n_z, n_p, alpha, color):
#    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,3), unpack=True)
#    dz = (max(z) - min(z)) / n_z
#    print (max(z), dz)
#    counter = 0
#    p_x = np.zeros(n_p)
#    I_y = np.zeros(n_p)
#    I_small_y = np.zeros(n_p)
#    for iz in range(n_z):
#        for ip in range(n_p):
#            p_x[ip] = p[counter]
#            I_y[ip] += abs(z[counter]) / f[counter]
#            if abs(z[counter]) > 2.0:
#                I_small_y[ip] += abs(z[counter]) / f[counter]
#            #print counter, z[counter]
#            counter = counter + 1
#    plt.plot(p_x, p_x**alpha * dz * I_y * 31., color=color)
#    plt.plot(p_x, p_x**alpha * dz * I_small_y * 31., color=color, linestyle=':')


#def read_jcr_at_z(filename, column, z_search, alpha, linestyle):
#    z, p, f = np.loadtxt(filename, skiprows=1, usecols=(0,1,column), unpack=True)
#    if (z_search < min(z) or z_search > max(z)):
#        print ("WARNING: z is out of range!")
#    best_z = -1000.0
#    for i in range(len(z)):
#        if abs(z[i] - z_search) < abs(best_z - z_search):
#            best_z = z[i]
#    ind = (z == best_z).nonzero()
#    p = p[ind] # GeV/c
#    E = p # GeV
#    f = f[ind] # m^-2 s^-1
#    EJ = f
#    print ("plot J in the range = ", min(EJ), max(EJ),np.mean(EJ))
#    plt.plot(E, E**alpha * EJ, linestyle)

#def analytical_solution_profile(H, D, vA): # kpc - 1e28 cm2/s - km/s
#    xi = 0.03 * vA * H / D
#    z = np.linspace(-H,H,100)
#    f = (1. - np.exp(-xi * (1. - abs(z) / H))) / (1. - np.exp(-xi))
#    plt.plot(z, f, 'k:')

#def kolmogorov_coefficient(D0, delta):
#    m_p = 0.938
#    p = np.logspace(-2,10,120)
#    E = np.sqrt(p * p + m_p * m_p)
#    beta = p / E
#    plt.plot(p, D0 * beta * p**delta, 'k:')

#def advection_timescale(vA, H):
#    p = np.logspace(-3,10,100)
#    t = ((H * 3.086e+19) / (vA * 1e3))  * (p / p)
#    plt.plot(p, t / 3.14e13, 'k:')

#def read_ratio_at_z(fn0, fn, column, z_search, alpha, linestyle):
#    z, p, f0 = np.loadtxt(fn0, skiprows=1, usecols=(0,1,column), unpack=True)
#    z, p, f  = np.loadtxt(fn,  skiprows=1, usecols=(0,1,column), unpack=True)
#    if (z_search < min(z) or z_search > max(z)):
#        print ("WARNING: z is out of range!")
#    best_z = -1000.0
#    for i in range(len(z)):
#        if abs(z[i] - z_search) < abs(best_z - z_search):
#            best_z = z[i]
#    ind = (z == best_z).nonzero()
#    p  = p[ind] # GeV/c
#    f0 = f0[ind] # GeV^-1 m^-2 s^-1
#    f  = f[ind]
#    print ("plot f in the range = ", min(f), max(f))
#    print ("plot f0 in the range = ", min(f0), max(f0))
#    plt.plot(p, (f - f0) / f0, linestyle)
