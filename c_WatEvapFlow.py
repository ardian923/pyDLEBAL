#! /usr/bin/env python
## Soil-Plant-Atmospheric Continuum calculation emphasizing on surface energy balance
## Developed initially by Ardiansyah (UNSOED), http://ardiansyah.net
## USE WITH "doubleLayerEB.xls" !! translated from EvPaddy-Laz9.26, CELIA STYLE
# write by Dr.Ardiansyah for research purpose
# ask for permission from me : ardi.plj@gmail.com (cc : ard@unsoed.ac.id) if you want to use this. It's free

import numpy as np  ##numerical python library
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sc  ##scientific python library
import math as mth  ##matemathical operation library

# we're going to use Microsoft Excel spreadsheet and xlrd xlwt module to link to Python
from xlrd import open_workbook
from xlwt import Workbook
from xlutils.copy import copy

# example writing and reading from/to excel
################################
workbook_name = '0_doubleLayerEB.xls'
sheet0 = 'Model'
sheet1 = 'Inputs'
sheet2 = '1-Layer'
sheet3 = '2-Layer'
sheet4 = 'Try'
sheet5 = 'Soil'
rb = open_workbook(workbook_name, formatting_info=True)
wb = copy(rb)

# for writing to excel
ws5 = wb.get_sheet(5)  # sheet indexing starts from 0
# ws5.write(0,0,'A1')
# wb.save(workbook_name)

# for reading from excel
rs1 = rb.sheet_by_name(sheet1)
# rs1.cell_value(0,0)  		# row and column index starts from 0
rs2 = rb.sheet_by_name(sheet2)


def get_cell_range(start_col, start_row, end_col, end_row):
    return [sheet.row_slice(row, start_colx=start_col, end_colx=end_col + 1) for row in xrange(start_row, end_row + 1)]

# BISMILLAHIRRAHMANIRRAHIIM
# Discretization


# READ DATA
depth = rs1.cell_value(115, 2)
bd = rs1.cell_value(116, 2)
n = rs1.cell_value(117, 2);
n = int(n)  ## number of nodes from 0 to 16
eps = rs1.cell_value(118, 2)
dz1 = rs1.cell_value(43, 2)
## Hydraulics parameters
ws = rs1.cell_value(121, 2)  ## saturated water content
ks = rs1.cell_value(122, 2)  ## saturated hydraulic conductivity, mm/s
ae = rs1.cell_value(123, 2)  ## air entry potential
a1 = rs1.cell_value(124, 2)
b1 = rs1.cell_value(125, 2)
c1 = rs1.cell_value(126, 2)
## Water retention parameters
a = rs1.cell_value(129, 2)
b = rs1.cell_value(130, 2)
c = rs1.cell_value(131, 2)
d = rs1.cell_value(132, 2)
e = rs1.cell_value(133, 2)
f = rs1.cell_value(134, 2)
my = rs1.cell_value(135, 2)  ## clay fraction of soil
## Simulation parameters
dt = rs1.cell_value(138, 2);
dt = float(dt)  ## timestep
inittime = rs1.cell_value(139, 2)
endtime = rs1.cell_value(140, 2)
im     = rs1.cell_value(141, 2)
maxits = rs1.cell_value(142, 2)
## Vapor related parameters
Mw     = rs1.cell_value(53, 2)  ## ava. in EB prog
R      = rs1.cell_value(48, 2)  ## ava. in EB prog
Dv     = rs1.cell_value(57, 2)  ## ava. in EB prog
rho_w  = rs1.cell_value(54, 2)  ## ava. in EB prog
p_init = rs1.cell_value(42, 2)  ## in kPa
#p_init = 1 * np.abs(p_init) * 0.102  ## convert kPa to mH2O
alt    = rs1.cell_value(12, 2)  ## ava. in EB prog
## Flux parameter for constant simulation
ETp  = rs1.cell_value(143, 2)  ## flux is converted in spreedsheet, from mm/day to m/(10 minutes), timestep is in 10 minutes
flux = rs1.cell_value(144, 2)  ## flux is converted in spreedsheet, from mm/day to m/(10 minutes), timestep is in 10 minutes
psurface = rs1.cell_value(145, 2)  ## p at surface =1, set upper boundary condition for evaporation (positive). psurface >=0 flux boundary condition

def geoGrid(depth, n, dz1):  ## Geometric grid with dry layer
    m = n - 2
    z = np.zeros(n + 1)
    sc = 0
    for i in xrange(1, m + 1):
        sc = sc + i * i
    dz = depth / sc
    z[0] = z[1] = 0
    if (dz1 != 0):
        z[2] = z[1] + dz1
        for i in xrange(2, m, 1):
            z[i + 1] = z[i] + dz * i * i
    else:
        for i in xrange(1, m, 1):
            z[i + 1] = z[i] + dz * i * i
    z[m + 1] = 1e+20  ##very big space = zero flux
    return z


def linGrid(depth, n, dz1):  ## Linear grid with dry layer
    m = n - 2
    z = np.zeros(n + 1)
    dz = depth / m
    z[0] = z[1] = 0
    if (dz1 != 0):
        z[2] = z[1] + dz1
        for i in xrange(2, m, 1):
            z[i + 1] = z[i] + dz
    else:
        for i in xrange(1, m, 1):
            z[i + 1] = z[i] + dz
    z[m + 1] = 1e+20  ##very big space = zero flux
    return z


## Liquid phase parameters to calculate
def waterCont(p, a, b, c, d, e, f):
    ## a, b, c, d, e and f is input in spreadsheet
    # p[p > ae] = ae #EVALUATE !!!!!!!!!!!!!!!!
    # define local variable pp for this function only
    #pp = np.zeros(p.size)
    pp = 1 * np.abs(p) * 0.102 # in mH2O, from kPa. Equation requires mH2O
    w = a / np.power((1 + np.power((b * pp), c)), d) + e * (1 - (np.log10(pp + 1) / f))
    u = np.power((1 + np.power((b * pp), c)), d)
    dwdp = (b * c * d * u * np.power((b * pp), c)) / (b * pp * np.power(u, 3)) - e * (1 / (f * (pp + 1)))
    return w, dwdp


def element_w():  ## wu in node i, wl in none i+1, average in the middle or in element
    w = 0.5 * (wnu + wnl)
    return w


def hydrConductivity(ks, wc, ws, a1, b1, c1):
    ##consult spreedsheet for a1, b1, and c1
    ##hydraulic conductivity is for layer or between node or element
    x1 = np.power((wc / ws), b1)
    x2 = 1 - x1
    if (wc == ws):
        kh = ks
    else:
        kh = ks * np.exp(-a1 * (np.power(x2, c1)))
    return kh


## ## Vapor phase parameters to calculate
def soilHumidity(p, T):
    #p still in kPa
    p = np.abs(p);
    p = -p  # required p to be (-) to obtain ha between 0 - 1
    soilHumidity = np.exp(Mw * p / (R * (T + 273)))  ##numpy array inserted
    return soilHumidity


def SVC(T):  ## Saturated Vapor Concentration, Kg/m3 T Kelvin
    Tk = T + 273
    rhov_sat = (1 / Tk) * np.exp(31.3716 - (6014.79 / Tk) - 0.00792495 * Tk)  ## g/m3
    rhov_sat = rhov_sat / 1000.  ## convert from g/m3 to Kg/m3
    drhov_sat = ((-1 / np.power(Tk, 2)) * rhov_sat) + (rhov_sat * ((6014.79 / np.power(Tk, 2)) - 0.00792495))
    drhov_sat = drhov_sat / 1000.  ## convert from g/m3 to Kg/m3
    return rhov_sat, drhov_sat


def AtmPressure(alt):  ## p_atm, in kPa, altitude in meter, ava. in EB prog
    AtmPressure = 101.3 * np.exp(-alt / 8200)
    return AtmPressure


def SVP(T):  ## Saturated Vapor Pressure, ava. in EB prog
    e_sat = 0.611 * np.exp(17.502 * T / (T + 240.97))  ## 3.8, Campbell 1998 (kPa)
    de_sat = 17.502 * 240.97 * e_sat / np.sqrt(240.97 + T)  ## 3.9, Campbell 1998
    return e_sat, de_sat  ##after return we repeat the name of variable/s so that we can assign it later to another new variable


def AVP(T, Twb, alt):  ## Actual Vapor Pressure, ava. in EB prog
    e_w, de_w = SVP(Twb)
    Pa = AtmPressure(alt) * 10  ## in hPa
    e_act = e_w - Pa * (T - Twb) * 0.00066 * (1 + (0.00115 * Twb))
    return e_act


def LatentHeatVaporize(T):  ##lambda, Latent heat of vaporization J/kg, ava. in EB prog
    LatentHeatVaporize = (45144 - 48 * T) / Mw
    return LatentHeatVaporize


### Liquid phase of water==================================
def k_bar(i):
    ##k average in time, water content calculated in node not in layer
    ## wat content in node, but in layer has subscript i+1/2
    eps = 0.5
    wi = 0.5 * (wiu[i] + wil[i])  ## upper at i itu = i, lower at i itu = i+1
    wn = 0.5 * (wnu[i] + wnl[i])
    ki = hydrConductivity(ks, wi, ws, a1, b1, c1)
    kn = hydrConductivity(ks, wn, ws, a1, b1, c1)
    kbar = (1 - eps) * ki + eps * kn  ## time averaged k, not used in this calculation
    return ki, kn


def q_liquid(i):
    ## i is element number; upper p is i, lower is i+1
    ki, kn = k_bar(i)
    qli = -(ki / (z[i + 1] - z[i]) * (pi[i + 1] - pi[i])) + ki   #tetap pada satuan kg/(m^2 s) ~ mm/s
    qln = -(kn / (z[i + 1] - z[i]) * (pn[i + 1] - pn[i])) + kn   #tetap pada satuan kg/(m^2 s) ~ mm/s
    J_liquid = (1 - eps) * qli + eps * qln  ##flux hasil perhitungan antara dua timestep, bukan dua level iterasi
    print "i, qli,      qln,      p[i+1],      p[i],      ki,      kn"
    print i, qli, qln, pi[i + 1], pi[i], ki, kn
    return qli, qln, J_liquid


def liqCoeff(i):
    wn = 0.5 * (wnu[i] + wnl[i])
    kn = hydrConductivity(ks, wn, ws, a1, b1, c1)
    UpLCoefn = eps * (kn / (z[i + 1] - z[i]))
    LowLCoefn = -eps * (kn / (z[i + 1] - z[i]))
    ResLCoefn = eps * kn
    return UpLCoefn, LowLCoefn, ResLCoefn


### Vapor phase of water==================================
def kv_kvT(i):
    c3 = 1 + 2.64 / np.sqrt(my);
    eta = 9.5 + 3 * w[i] - (9.5 - 1) * np.exp(-np.power((c3 * w[i]), 4)) #conductivity (W m-1 K-1), eq.4.20, Campbell
    ## humidity initial, and final end of timestep
    # value of p is already in kPa
    h[i] = soilHumidity(p[i], T[i])       #no unit
    hi[i] = soilHumidity(pi[i], Ti[i])    #no unit
    hn[i] = soilHumidity(pn[i], Tn[i])    #no unit
    hbar_i = 0.5 * (hi[i] + hi[i + 1])
    hbar_n = 0.5 * (hn[i] + hn[i + 1])
    ## porosity
    phi_i = ws - 0.5 * (wnu[i] + wnl[i])  #no unit
    phi_n = ws - 0.5 * (wiu[i] + wil[i])  #no unit
    ## vapor concentration
    rhov_sati, drhov_sati = SVC(Ti[i]);   #Kg/m3
    rhov_i = h[i] * rhov_sati  ## actual vapor concentration at i (iteration level)
    rhov_satn, drhov_satn = SVC(Tn[i]);
    rhov_n = h[i] * rhov_satn  ## actual vapor concentration at n (iteration level)
    ## conductivity for vapor, isothermal and non-isothermal(tempereture influenced) conductivity
    kvi = 0.66 * phi_i * Dv * (Mw / (R * (Ti[i] + 273))) * rhov_i  ## Dv, Mw inputted in spreedsheet
    kvTi = 0.66 * phi_i * Dv * drhov_sati * eta * hbar_i
    kvn = 0.66 * phi_n * Dv * (Mw / (R * (Tn[i] + 273))) * rhov_n
    kvTn = 0.66 * phi_n * Dv * drhov_satn * eta * hbar_n
    ## convert to m/s
    #kvi = kvi * (1 / 0.102) * (1 / rho_w)  ## 1 kPa = 0.102 mH2O, convert matriks to kPa before using this conduct
    #kvTi = kvTi * (1 / 0.102) * (1 / rho_w)  ## 1 kPa = 0.102 mH2O, convert matriks to kPa before using this conduct
    #kvn = kvn * (1 / 0.102) * (1 / rho_w)  ## 1 kPa = 0.102 mH2O, convert matriks to kPa before using this conduct
    #kvTn = kvTn * (1 / 0.102) * (1 / rho_w)  ## 1 kPa = 0.102 mH2O, convert matriks to kPa before using this conduct
    return kvi, kvn, kvTi, kvTn


def q_vapor(i):  ## non-isothermal flow, in case of isothermal deltaT = 0
    lamda = LatentHeatVaporize(T)  ## J/kg
    kvi, kvn, kvTi, kvTn = kv_kvT(i)
    qvi = (1 / (z[i + 1] - z[i])) * (-kvi * (pi[i + 1] - pi[i]) - kvTi * (Ti[i + 1] - Ti[i]))
    qvn = (1 / (z[i + 1] - z[i])) * (-kvn * (pn[i + 1] - pn[i]) - kvTn * (Tn[i + 1] - Tn[i]))
    J_vapor = (1 - eps) * qvi + eps * qvn  ## flux hasil perhitungan antara dua timestep, bukan dua level iterasi
    ## htsrc = 10 * lamda * (kvn * rho_w) * (pn[i+1]-pn[i]) ## input for heat flux calculation(J/m2.s), qT = q + htsrc
    print "i,      qvi,       qvn,      kvi,      kvn,      kvTi,      kvTn"
    print i, qvi, qvn, kvi, kvn, kvTi, kvTn
    return qvi, qvn, J_vapor


def vapCoeff(i):
    kvi, kvn, kvTi, kvTn = kv_kvT(i)
    UpVCoefn = eps * (kvn / (z[i + 1] - z[i]))
    LowVCoefn = -eps * (kvn / (z[i + 1] - z[i]))
    ResVCoefn = -eps * (kvTn * (Tn[i + 1] - Tn[i])) / (z[i + 1] - z[i])
    return UpVCoefn, LowVCoefn, ResVCoefn


### Root extraction from soil ==================================
### Input root extraction in soil is also transpiration rate

def initSoilCondition(T):
    m = n - 2  ## 15 end node number that important
    p = np.zeros(n)
    p[:] = p_init
    pi = np.zeros(n)
    pi[:] = p_init
    pn = np.zeros(n)
    pn[:] = p_init
    ## ks, ws, pe are global variables that includ
    wu, dwudp = waterCont(p, a, b, c, d, e, f);
    wl = 1 * wu
    dwnudp = 1 * dwudp
    wiu = 1 * wu
    wil = 1 * wu
    wnu = 1 * wu
    wnl = 1 * wu
    wnu[m + 1] = wnl[m + 1] = wu[m]
    # p[0] = p[1]
    # pi = pn = p ## pi initial timestep, p middle of timestep, pn end of timestep
    h = soilHumidity(p, T);     #no unit
    hi = soilHumidity(pi, T);   #no unit
    hn = soilHumidity(pn, T);   #no unit
    z[0] = 0  # -1E+20;
    z[m + 1] = +1E+20  ## No upward vapor flux into bottom
    for i in xrange(1, m + 2, 1):
        v[i] = (z[i + 1] - z[i - 1]) / 2  ## dz
    ## Tn[:] = 26 ## for isothermal soil, for non isothermal coupled with soil heat flow
    ## Ti = 1 * Tn ##remove when this program coupled with Soil heat flow (non-isothermal)
    return wu, wl, wiu, wil, wnu, wnl, dwudp, dwnudp, pi, p, pn, hi, h, hn, v


def boundaryCondition(flux, evap, psurface):
    Jl = np.zeros(n)
    qli = np.zeros(n)
    qln = np.zeros(n)  ## assign liquid flux in all nodes
    Jv = np.zeros(n)
    qvi = np.zeros(n)
    qvn = np.zeros(n)  ## assign vapor flux in all nodes
    UpLCoefn = np.zeros(n)
    LowLCoefn = np.zeros(n)
    ResLCoefn = np.zeros(n)
    UpVCoefn = np.zeros(n)
    LowVCoefn = np.zeros(n)
    ResVCoefn = np.zeros(n)
    ## surface boundary condition
    if (psurface < 0):
        if (psurface > ae):
            p[1] = ae
        else:
            p[1] = psurface
    else:
        ## flux boundary condition
        Jl[0] = flux
        qli[0] = qln[0] = Jl[0]  ## flux boundary condition, time average flux for flux boundary condition
        Jv[0] = evap
        qvi[0] = qvn[0] = Jv[0]
        UpLCoefn[0] = LowLCoefn[0] = ResLCoefn[0] = 0
        UpVCoefn[0] = LowVCoefn[0] = ResVCoefn[0] = 0
    return Jl, qli, qln, Jv, qvi, qvn, UpLCoefn, LowLCoefn, ResLCoefn, UpVCoefn, LowVCoefn, ResVCoefn


def thomasAlgorithm(i1, A, B, C, D):  ##tridiagonal matrix solution
    m = n - 2  ## 15 nodes to evaluate, n = 17, including 0 and 17th
    x = np.zeros(n)
    for i in xrange(i1, m):
        C[i] = C[i] / B[i]  ## update C
        D[i] = D[i] / B[i]  ## update D
        B[i + 1] = B[i + 1] - A[i + 1] * C[i]  ## update B
        D[i + 1] = D[i + 1] - A[i + 1] * D[i]  ## update D
        #print "Thomas A B C D", A[i], B[i], C[i], D[i]
    x[m] = D[m] / B[m];
    #print "PPPPPPPPPPPPPPPPPPPPPPPPPPP", pn[m], D[m] / B[m]
    #print x[15]
    for i in xrange(m - 1, i1 - 1, -1):
        x[i] = D[i] - C[i] * x[i + 1]
    return x


def solverWatEvapFlow(n, dt, flux, evap, psurface):
    ## use below's global variable
    global wu, wl, wiu, wil, wnu, wnl, dwudp, dwnudp, pi, p, pn, hi, h, hn, v
    global se, nits
    ## set local variables for computing
    cpu = np.zeros(n)
    m = n - 2  ## 15 nodes to evaluate, n = 17, including 0 and 17th
    A = np.zeros(m + 1);
    B = np.zeros(m + 1);
    C = np.zeros(m + 1);
    D = np.zeros(m + 1)  ## tridiagonal matrix coefficient
    ## Apply boundary condition while define local variables
    Jl, qli, qln, Jv, qvi, qvn, UpLCoefn, LowLCoefn, ResLCoefn, UpVCoefn, LowVCoefn, ResVCoefn = boundaryCondition(flux,
                                                                                                                   evap,
                                                                                                                   psurface)
    se = 1;
    nits = 0  ## number of iteration until solution convergence
    while (se > im):  ## loop until mass balance error less than tolerance
        nits = nits + 1  ## update number of iteration
        if nits >= maxits:
            break  ## stop iterating on inconvergence
        p = 1 * pn  ## assign previous iteration to current intial conditon, untuk evaluasi error
        for i in xrange(1, m + 1):  ## i = 1 to 15
            cpu[i] = v[i] * dwnudp[i] / dt  ## Water capacity
            ## Liquid flux and liquid coefficient
            qli[i], qln[i], Jl[i] = q_liquid(i)
            UpLCoefn[i], LowLCoefn[i], ResLCoefn[i] = liqCoeff(i)
            ## Vapor flux and vapor coefficient
            qvi[i], qvn[i], Jv[i] = q_vapor(i)
            UpVCoefn[i], LowVCoefn[i], ResVCoefn[i] = vapCoeff(i)
            ## obtain tridiagonal matrix coefficient
            A[i] = UpLCoefn[i - 1] + UpVCoefn[i - 1]
            C[i] = -LowLCoefn[i] - LowVCoefn[i]
            B[i] = -UpLCoefn[i] + LowLCoefn[i - 1] - UpVCoefn[i] + LowVCoefn[i - 1] - cpu[i]
            if (i == 1):
                D[i] = -(1 - eps) * ((qli[i - 1] - qli[i]) + (qvi[i - 1] + qvi[i])) - cpu[i] * p[i] + (
                            wu[i] - wiu[i]) * (v[i] / dt) - (ResLCoefn[i] - ResLCoefn[i - 1]) - \
                       (ResVCoefn[i] - ResVCoefn[i - 1]) - eps * (qln[i - 1] + qvn[i - 1])
            else:
                D[i] = -(1 - eps) * ((qli[i - 1] - qli[i]) + (qvi[i - 1] + qvi[i])) - cpu[i] * p[i] + (
                            wu[i] - wiu[i]) * (v[i] / dt) - (ResLCoefn[i] - ResLCoefn[i - 1]) - (
                                   ResVCoefn[i] + ResVCoefn[i - 1])
        ## preparing to solve tridiagonal matrix
        if (psurface < 0):  ## dirichlet BC not flux BC, start calculation from i = 2
            D[1] = 0;
            C[1] = 0;
            i1 = 2
        else:
            i1 = 1
        ## solve tridiagonal matrix
        pn = 1 * thomasAlgorithm(i1, A, B, C, D)
        ## calculate water content at each node and total error of mass balance
        wu, dwudp = waterCont(p, a, b, c, d, e, f)  ## if p is an array, wu and dwudp are also array, no need looping
        wnu, dwnudp = waterCont(pn, a, b, c, d, e, f)
        #print "pi ", pi
        print "p ", p
        print "pn", pn
        print "wu", wu
        print "wnu", wnu
        print "Jv[0]", Jv[0]
        print "Jl[0]", Jl[0]
        se = 0
        for i in xrange(1, m + 1):  ## i = 1 to 15
            wl[i] = 1 * wu[i + 1];
            wnl[i] = 1 * wnu[i + 1];
            se = se + abs(wu[i] - wnu[i]);
            print "ERROR EUY !! = ", se
        ## while loop ends here
    ## when convergence, calculate water depth and assign end of timestep variable to intial for next initial iteration
    pn[m + 1] = pn[m]  ## Unit gradient drainage at bottom, no drainage
    if (nits <= maxits):
        sw = 0
        for i in xrange(m, 1, -1): sw = sw + v[i] * (
                    wnu[i] - wiu[i] + wnl[i - 1] - wil[i - 1]) / 2  ## water depth storage in soil depth
        ## assign result of wnu to wu and wiu
        wu = 1 * wnl;
        wiu = 1 * wnl
        wl = 1 * wnl;
        wil = 1 * wnl
        pi = 1 * pn  ## memperkenalkan "pi" untuk assign pn ke p timestep sebelumnya yang bisa digunakan dalam menghitung flux
        success = True
    else:
        success = False
    w = element_w()
    return se, nits  ## end of solver here


####### MAIN PROGRAM #######
if __name__ == '__main__':  ## Run as standalone program
    ## Define numerical variables
    m = n - 2
    # z = geoGrid(depth, n , dz1) ## geometric grid
    z = linGrid(depth, n, dz1)  ## linear grid
    for col in xrange(1, m + 2, 1):
        ws5.write(5, 7 + col, z[col])

    v = np.zeros(n)  ## dz
    ## htsrc = np.zeros(n) ## heat source for heat flow problem
    # se = 0 ## initial starter to run "while" loop in solverWatEvapFlow
    # nits = 0 ## sum of error and number of iteration to convergence
    if (ae > 0): ae = -ae  ## potential is (-) : suction

    T = np.zeros(n);
    T[:] = 26;
    Ti = 1 * T;
    Tn = 1 * T  ## When soil is non-isothermal, remove this line to replace T with temperature profile
    wu, wl, wiu, wil, wnu, wnl, dwudp, dwnudp, pi, p, pn, hi, h, hn, v = initSoilCondition(T)
    w = element_w()

    dt = dt / 3600  ##

    time = np.arange(inittime, endtime + 0.1, dt)  ##create time array from 0 to endtime
    ## ==================== Atmospheric Factor Influencing Evaporation============
    # obtain data value of Ta and Twb for calculating RH of air  # ava. in EB prog
    ## sensor measurement input or Fourier input or data input # ava. in EB prog
    Ta = np.zeros(time.size)  ##allocating space for Ta array # ava. in EB prog
    Twb = np.zeros(time.size)  ##allocating space for Ta wet bulb array # ava. in EB prog
    RH = np.zeros(time.size)  ##allocating space for RH array, calculated from Ta wet bulb # ava. in EB prog
    i = 0
    for row in range( 12, int(12+1 +(6*endtime)) ): # endtime = 24 or 48 hours or more
        Ta[i] = rs2.cell_value(row, 7)
        Twb[i] = rs2.cell_value(row, 8)
        i = i + 1

    ## calculate air humidity from data # ava. in EB prog
    e_sat, de_sat = SVP(Ta)
    e_act = AVP(Ta, Twb, alt) #PERHATIKAN BARIS INI, RH AKAN BERNILAI 1, MENYEBABKAN PERHITUNGAN EVAP ADA PEMBAGI NOL
    RH = e_act / e_sat
    print 'RH', RH

    ## solve for water-evap flow every timestep
    for i in xrange(time.size):
        print "=========================================================================="
        print "time = ", time[i], "======================================================"
        print "Soil depth discretization", z

        ## for constant evaporation and infiltration through day, single value inputted. For not changing value put in spreedsheet and
        #evap = (ETp / (1 - RH[i])) * (h[1] - RH[i])  ## actual evaporation from bare surface
        evap = ETp #try, evaporation were calculated separately, not from humidity difference between soil and air
        ## obtain the value through looping
        se, nits = solverWatEvapFlow(n, dt, flux, evap, psurface)
        # write simulation result to spreedsheet
        ws5.write(7 + i, 4, time[i])
        ws5.write(7 + i, 5, nits)
        ws5.write(7 + i, 6, se)
        ws5.write(7 + i, 7, evap)
        ws5.write(7 + i, 8, flux)
        ##    for row in time:
        for col in xrange(1, m + 1, 1):
            ws5.write(7 + i, col + 8, wnu[col])
            # ws5.write(7+i,6, col+8, hn[col])
        ## dynamic graph in matplotlib --> not implemented
        # static graph for selected time
        if (int(time[i]) == 0) and (time[i]-int(time[i])<0.1):
            string_time0 = 'time = 0'
            z0 = z[1:18]       # slicing z vector to make same size as p and w
            p0 = pn
            w0 = w
        if (int(time[i]) == 2) and (time[i]-int(time[i])<0.1):
            string_time2 = 'time = 2'
            p2 = pn
            w2 = w
        if (int(time[i]) == 10) and (time[i]-int(time[i])<0.1):
            string_time10 = 'time = 10'
            p10 = pn
            w10 = w
        if (int(time[i]) == 24) and (time[i]-int(time[i])<0.1):
            string_time24 = 'time = 24'
            p24 = pn
            w24 = w

    # Save editted worksheet
    wb.save(workbook_name)

    # plot results
    mpl.rc('font', family='serif')
    mpl.rc('font', serif='Helvetica Neue')
    mpl.rc('text', usetex='true')
    mpl.rcParams.update({'font.size': 16})

    fig1, (ax1, ax2) = plt.subplots(ncols=2)
    ax1.plot(p0, z0, 'r-', label=string_time0)
    ax1.plot(p2, z0, 'k-', label=string_time2)
    ax1.plot(p10,z0, 'g-', label=string_time10)
    ax1.plot(p24,z0, 'b-', label=string_time24)
    ax1.legend(loc='upper left', shadow=False, frameon=False)
    ax1.invert_yaxis()
    ax1.set_xlim([-1500, 0])
    #plt.ylim(ymax = 1, ymin = 0)
    #plt.title('Matric Suction')
    ax1.set_xlabel('Matric Potential ($\psi - kPa$)')
    ax1.set_ylabel('Soil Depth ($m$)')

    ax2.plot(w0, z0, 'r-', label=string_time0)
    ax2.plot(w2, z0, 'k-', label=string_time2)
    ax2.plot(w10,z0, 'g-', label=string_time10)
    ax2.plot(w24,z0, 'b-', label=string_time24)
    ax2.legend(loc='upper left', shadow=False, frameon=False)
    ax2.invert_yaxis()
    ax2.set_xlim([0, 0.45])
    #plt.ylim(ymax = 1, ymin = 0)
    #plt.title('Matric Suction')
    ax2.set_xlabel('Water Content ($\\theta - cm^3/cm^3$)')
    ax2.set_ylabel('Soil Depth ($m$)')

    plt.show()
