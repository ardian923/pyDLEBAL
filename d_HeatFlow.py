#! /usr/bin/env python
## Soil-Plant-Atmospheric Continuum calculation emphasizing on surface energy balance
## Developed initially by Ardiansyah (UNSOED), http://ardiansyah.net
## USE WITH "OoPython-ETo.ods" !! 

##ardiansyah@AL-FATIH-II:~/Desktop/OoPython-ETo$ oocalc OoPython-ETo.ods -accept="socket,host=localhost,port=2002;urp;StarOffice.ServiceManager"
##python

import os as os ##operating system function
from oosheet import OOSheet as pyo ##open office sheet
import numpy as np ##numerical python library
import pylab as plt ##plotting library
import scipy as sc ##scientific python library
import math as mth ##matemathical operation library

import FuncLibrary as fl

## BISMILLAHIRRAHMANIRRAHIIM

## Discretization
depth = pyo('Inputs.c116').value
bd    = pyo('Inputs.c117').value  ## bulk density
n     = pyo('Inputs.c118').value; n = int(n) ## number of nodes from 0 to 16
eps   = pyo('Inputs.c119').value
dz1   = pyo('Inputs.c44').value ## Dry layer depth
wc    = ????? ##water content 
## 

def geoGrid(depth, n, dz1): ## Geometric grid with dry layer
    m = n - 2
    z = np.zeros(n+1)
    sc = 0
    for i in xrange(1,m+1):
        sc = sc + i*i
    dz = depth/sc
    z[0] = z[1] = 0
    if (dz1<>0):
        z[2] = z[1] + dz1
        for i in xrange(2,m,1):
            z[i+1] = z[i] + dz*i*i
    else:
        for i in xrange(1,m,1):
            z[i+1] = z[i] + dz*i*i
    z[m+1] = 1e+20 ##very big space = zero flux
    return z

def linGrid(depth, n, dz1): ## Linear grid with dry layer
    m = n - 2
    z = np.zeros(n+1)
    dz = depth/m
    z[0] = z[1] = 0
    if (dz1<>0):
        z[2] = z[1] + dz1
        for i in xrange(2,m,1):
            z[i+1] = z[i] + dz
    else: 
        for i in xrange(1,m,1):
            z[i+1] = z[i] + dz 
    z[m+1] = 1e+20 ##very big space = zero flux
    return z
 
## Soil thermal properties : heat capacity and thermal conductivity
def heatCapacity(bd, wc):
    heatCap = 2400000*bd/2.65 + 4180000*wc
    return heatCap

def thermalConductivity(i, bd,wc,my,Ti):
    c1 = 0.65 - 0.78*bd + 0.6*sqr(bd)
    c2 = 1.07*bd
    c3 = 1+2.64/sqrt(my) ## my = clay fraction
    c4 = 0.03+0.1*sqr(bd)
    eta= 9.5+6*wc-(9.5-1)*exp(-power((c3*wc),4))
    porosity = ws[i] - 0.5*(wnu[i]+wnl[i])
    hr = Humidity(i,p[i])
    s  = dSatVap(Ti)/1000
    L  = LatentHeatVaporize(Ts)/(Mw*1000) ## J/g
    ThermalConductivity = c1+c2*wc-(c1-c4)*exp(-power((c3*wc),4))+0.66*Dv*porosity*hr*s*eta*L
    return ThermalConductivity

## Initial condition
def initSoilTemp(Ta, Tb):
    m = n-2 ## 15 end node number that important
    Tn[0] = Ta
    for i in xrange(1,m+2,1):
        Tn[i] = Tb
    return Tn

def updateSoilTemp(t, m)
    Tn[0] = Ta
    Tn[1] = Ts
    for i in xrange(0, m+2, 1):
        Ti[i] = Tn[i]
    return Tn, Ti

def q_heat(i): 
    Tbar[i]   = (1-eps)*Ti[i] + eps*Tn[i]
    Tbar[i+1] = (1-eps)*Ti[i+1] + eps*Tn[i+1]

    J_heati   = -eps*kT[i]/(z[i+1]-z[i])*(Ti[i+1]-Ti[i])
    J_heat    = -kT[i]/(z[i+1]-z[i])*(Tbar[i+1]-Tbar[i])
    dJhdTu    = eps*kT[i]/(z[i+1]-z[i])
    dJhdTl    = eps* -kT[i]/(z[i+1]-z[i])
    return J_heati, J_heat, dJhdTu, dJhdTl

def thomasAlgorithm(i1, A, B, C, D): ##tridiagonal matrix solution
    m = n-2 ## 15 nodes to evaluate, n = 17, including 0 and 17th
    for i in xrange(i1, m):
        C[i] = C[i]/B[i] ## update C 
        D[i] = D[i]/B[i] ## update D
        B[i+1] = B[i+1] - A[i+1]*C[i] ## update B
        D[i+1] = D[i+1] - A[i+1]*D[i] ## update D

    D[m] = D[m]/B[m]
    for i in xrange(m-1, i1-1, -1):
        D[i] = D[i] - C[i] * D[i+1]
    return D

def q_vapor(i): ## non-isothermal flow, in case of isothermal deltaT = 0
    lamda = LatentHeatVaporize(T) ## J/kg
    kvi, kvn, kvTi, kvTn = kv_kvT(i)
    qvi = (1/(z[i+1]-z[i]))*(-kvi*(pi[i+1]-pi[i]) - kvTi*(Ti[i+1]-Ti[i]))
    qvn = (1/(z[i+1]-z[i]))*(-kvn*(pn[i+1]-pn[i]) - kvTn*(Tn[i+1]-Tn[i]))
    J_vapor = (1-eps)*qvi + eps*qvn ## flux hasil perhitungan antara dua timestep, bukan dua level iterasi
    ## htsrc = 10 * lamda * (kvn * rho_w) * (pn[i+1]-pn[i]) ## input for heat flux calculation(J/m2.s), qT = q + htsrc
    print "qvi, qvn", qvi, qvn
    return qvi, qvn, J_vapor

def solverHeatFlow(n, dt, bd, my, blcond, nbd, w):
    global Ti, Tn
    m = n-2 ## 15 nodes to evaluate, n = 17, including 0 and 17th
    for i in xrange(1, m+1):
        cp[i] = FWG_HC(wg[i])* 1e6 ## J/m3.K  ## ??????????????????????????????????????????????????????
        kT[i] = FWG_K(wg[i])                  ## ??????????????????????????????????????????????????????
        J_heati[i], J_heat[i], dJhdTu[i], dJhdTl[i] = q_heat(i)
        CC[i]   = dJhdTl[i]
        AA[i+1] = C[i];
        BB[i]   = -(dJhdTl[i-1]-dJhdTu[i]) + cp[i]*(z[i+1]-z[i-1])/(2*dt);
        DD[i]   = Jhi[i-1]-Jhi[i] + Ti[i]*cp[i]*(z[i+1]-z[i-1])/(2*dt) + htsrc[i];

    DD[1] = DD[1] ## + kT[0] * f*Tn[0] - Rn + LE;
    DD[m] = DD[m] + kT[m] * eps*Tn[m+1]
    de = DD[1]: beex = BB[1]: ce = CC[1]
    di = DD[2]: biix = BB[2]: ci = CC[2]: ai = AA[2]
    LEf = (di-(ai*Ti[1]+biix*Ti[2]+ci*Ti[3])) ## Dry layer evaporation if dz<>0
    if (dz1 = 0):
        DD[2] = DD[2]
    else:
        DD[2] = DD[2] + LEf
    
    i1 = 2 ## not change Tn[1] as surface boundary condition
    Tn = thomasAlgorithm(i1, AA, BB, CC, DD)
    global Ghf
    if (dz1 = 0):
        Ghf = (cp[1]*(z[2]-z[1])*(Tn[1]-Ti[1]))/dt + Jh[1]
    else:
        Ghf = (cp[1]*(z[2]-z[1])*(Tn[1]-Ti[1]))/dt + Jh[1] - (LEf)
    dGhfdTs = -(cp[1]*(z[2]-z[1]) / dt
    return Tn, Ghf

###############################################################
if __name__ == '__main__': # run heatflow simulation as standalone program
    ## Main Program
    ## Simulation parameters
    dt = pyo('Inputs.c139').value; dt = float(dt) ## timestep
    inittime = pyo('Inputs.c140').value
    endtime = pyo('Inputs.c141').value

    ## Define numerical variables    
    m = n - 2
    z = linGrid(depth, n , dz1) ## linear grid
    for col in xrange(1, m+2, 1):
        pyo('Soil.aa6')[0][col].value = z[col]

    # read input variables from spreadsheet
    Tb =  ## Initial soil temperature for each layer 
    Tn = initSoilTemp(Ta, Tb) ## Soil temperature profile 
    AA = np.zeros(n); BB = np.zeros(n); CC = np.zeros(n); DD = np.zeros(n) ## tridiagonal matrix coefficient
    
    ## solve for heat flow every timestep
    time = 0
    i = 0
    while (time <= endtime):
        Ta =  ## air temperature, measured value or fourier as function of time
        Ts =  ## surface temperature, CALCULATED, measured value or fourier as function of time
        RH =  ## RH, measured value or fourier as function of time
        Tn, Ti = updateSoilTemp(time, m)
        print "time = ", time, "======================================================"
        # write simulation result to spreedsheet
        pyo('Soil.g8')[i]['E'].value = time
        pyo('Soil.g8')[i]['I'].value = flux
        ## for row in time:
        for col in xrange(0, m+2, 1):
            pyo('Soil.aa8')[i][col].value = Tn[col]

        ## dynamic graph in matplotlib

        ## solve and obtain Tn value within time looping
        solverHeatFlow(n, dt, flux, evap, psurface) ## resulting new Tn value
    
        # update time step
        time = time + dt/3600.
        i = i + 1


