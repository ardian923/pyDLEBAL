import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#p = input('enter water potential :')

a = 0.343
b = 0.04133
c = 5.603
d = 0.8215
e = 0.08624
f = 14.37

def waterCont(p, a, b, c, d, e, f):
    ##consult spreedsheet for a, b, c, d, e and f
    # p[p > ae] = ae #EVALUATE !!!!!!!!!!!!!!!!
    #pp = np.zeros(p.size)
    pp = 1 * np.abs(p)  # in mH2O
    w = a / np.power((1 + np.power((b * pp), c)), d) + e * (1 - (np.log10(pp + 1) / f))
    u = np.power((1 + np.power((b * pp), c)), d)
    dwdp = (b * c * d * u * np.power((b * pp), c)) / (b * pp * np.power(u, 3)) - e * (1 / (f * (pp + 1)))
    return w, dwdp


# Hydraulic conductivity
a1 = 40.3
b1 = 0.672
c1 = 3.31
ks = 0.1093 # mm/s
ws = 0.43

def hydrConductivity(ks, wc, ws, a1, b1, c1):
    ##consult spreedsheet for a1, b1, and c1
    ##hydraulic conductivity is for layer or between node or element
    x1 = np.power((wc / ws), b1)
    x2 = 1 - x1
    kh = ks * np.exp(-a1 * (np.power(x2, c1)))
    #if (wc == ws):
    #    kh = ks
    #else:
    #    kh = ks * np.exp(-a1 * (np.power(x2, c1)))
    return kh

if __name__ == '__main__':  ## Run simulation as standalone program
    ## Main Program

    p_kPa = np.array([1, 10, 33, 100, 500, 1000, 1500])
    p_mmH2O = np.array([102.0, 1019.7, 3365.0, 10197.0, 50985.0, 101970.0, 152955.0])
    p_cmH2O = np.array([10.20, 101.97, 336.50, 1019.70, 5098.50, 10197.00, 15295.50])
    p_mH2O = np.array([0.1020, 1.0197, 3.3650, 10.1970, 50.9850, 101.9700, 152.9550])

    p = p_mH2O
    theta, dtheta = waterCont(p, a, b, c, d, e, f)
    wc = theta
    kh = hydrConductivity(ks, wc, ws, a1, b1, c1)
    print 'Theta = ', theta
    print 'dTheta= ', dtheta
    print 'Kh = ', kh

    # plot results
    mpl.rc('font', family='serif') 
    mpl.rc('font', serif='Helvetica Neue')
    mpl.rc('text', usetex='true')
    mpl.rcParams.update({'font.size': 16})

    # 1-Layer Result
    plt.figure(1)
    plt.subplot(211)
    plt.plot(theta, p, 'r-', label='Water Retention')
    plt.xlabel('$theta$ ($cm^3/cm^3$))')
    plt.ylabel('$\psi$ ($kPa$)')
    plt.legend(loc='upper left', shadow=False, frameon=False)

    plt.subplot(212)
    plt.plot(theta, kh, 'g-', label='Hydr Conductivity')
    plt.xlabel('$theta$ ($cm^3/cm^3$))')
    plt.ylabel('Kh ($mm/s$)')
    plt.legend(loc='upper left', shadow=False, frameon=False)

    plt.show()

