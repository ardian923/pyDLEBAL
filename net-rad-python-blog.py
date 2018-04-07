<pre lang="PYTHON" line="1">
#! /usr/bin/env python

import os as os ##operating system function
from oosheet import OOSheet as pyo ##open office sheet
import numpy as np ##numerical python library
import pylab as plt ##plotting library
import scipy as sc ##scientific python library
import math as mth ##matemathical operation library

##obtain all input from spreedsheet
##consult spreedsheet for parameters note, symbol, unit, etc
##Boundary Layer Parameters
LAI = pyo('Inputs.c8').value

##Radiation Parameters
alb_s = pyo('Inputs.c22').value
alb_c = pyo('Inputs.c23').value
eps_s = pyo('Inputs.c24').value
eps_c = pyo('Inputs.c25').value
deg = pyo('Inputs.c27').value
min = pyo('Inputs.c28').value
sec = pyo('Inputs.c29').value
hem = pyo('Inputs.c30').string
rad = pyo('Inputs.c31').value
daystrt = pyo('Inputs.c32').value

##Constant for fourier Ta
deltaTa = pyo('Inputs.c60').value
zTa = pyo('Inputs.c61').value
aTa0 = pyo('Inputs.c62').value
aTa = np.array([pyo('Inputs.c63').value, pyo('Inputs.c64').value, pyo('Inputs.c65').value, pyo('Inputs.c66').value])
bTa = np.array([pyo('Inputs.c67').value, pyo('Inputs.c68').value, pyo('Inputs.c69').value, pyo('Inputs.c70').value])

##Constant for fourier Tc
deltaTc = pyo('Inputs.c72').value
zTc = pyo('Inputs.c73').value
aTc0 = pyo('Inputs.c74').value
aTc = np.array([pyo('Inputs.c75').value, pyo('Inputs.c76').value, pyo('Inputs.c77').value, pyo('Inputs.c78').value])
bTc = np.array([pyo('Inputs.c79').value, pyo('Inputs.c80').value, pyo('Inputs.c81').value, pyo('Inputs.c82').value])

##Constant for fourier Ts
deltaTs = pyo('Inputs.c84').value
zTs = pyo('Inputs.c85').value
aTs0 = pyo('Inputs.c86').value
aTs = np.array([pyo('Inputs.c87').value, pyo('Inputs.c88').value, pyo('Inputs.c89').value, pyo('Inputs.c90').value])
bTs = np.array([pyo('Inputs.c91').value, pyo('Inputs.c92').value, pyo('Inputs.c93').value, pyo('Inputs.c94').value])

##Constant for fourier RH
deltaRH = pyo('Inputs.c96').value
zRH = pyo('Inputs.c97').value
aRH0 = pyo('Inputs.c98').value
aRH = np.array([pyo('Inputs.c99').value, pyo('Inputs.c100').value, pyo('Inputs.c101').value, pyo('Inputs.c102').value])
bRH = np.array([pyo('Inputs.c103').value, pyo('Inputs.c104').value, pyo('Inputs.c105').value, pyo('Inputs.c106').value])

sb = pyo('Inputs.c56').value

os.system('clear') #clear screen on linux

def SVP(T): ## Saturated Vapor Pressure, T in celcius
        e_sat  = 0.611*np.exp(17.502*T/(T+240.97)) ## 3.8, Campbell 1998 (kPa)
        de_sat = 17.502*240.97*e_sat/np.sqrt(240.97+T)   ## 3.9, Campbell 1998
        return e_sat, de_sat ##after return we repeat the name of variable/s so that we can assign it later to another new variable

def fourier(delta, j, z, a0, a, b): ##create fourier function for air temperature and air humidity
    i = np.arange(1, 5) ##create i array for sigma 
    w = 2 * np.pi / 24
    sigma = np.exp(-z*(i*w/(2*delta))**0.5) * (a * np.cos(i*w*time[j] - z*(i*w/(2*delta))**0.5)) + (b * np.sin(i*w*time[j] - z*(i*w/(2*delta))**0.5))
    fourier = a0 + sum(sigma)
    return fourier

def ConvertLattitude(deg, min, sec, hem):
    LattRadian = deg+(min/60.0)+(sec/3600.0)
    LattRadian = (np.pi/180)*LattRadian
    if hem == 'south':
        LattRadian =-LattRadian 
    else:
        LattRadian = LattRadian
    return LattRadian

###########################################################################
def Radiation(daystrt, hour, rad, Ta, RH): ## Calculate Radiation using input hour
    ## rad is total daily measured global radiation (MJ/m2)
    la = ConvertLattitude(deg, min, sec, hem)
    sin_la = np.sin(la) ## sin lattitude
    cos_la = np.cos(la) ## cos lattitude
    sin_d = 0.3985*np.sin(4.869+0.0172*daystrt+0.03345*np.sin(6.224+0.0172*daystrt))  ## sin declination 12.6, Campbell
    cos_d = np.sqrt(1-np.power(sin_d, 2))
    ## Prepare to calculate Total Radiation
    hrs = np.arange(1., 25., 1)
    sin_e = np.zeros(hrs.size)
    sin_e = sin_d*sin_la + cos_d*cos_la*np.cos(0.2618*(hrs-12))  ## Solar elevation angle 12.8, returns array of sin_e
    RadTot = 0
    for i in xrange(hrs.size):
        if sin_e[i]>0:
            RadTot = RadTot + 1360*sin_e[i] ## W-hr/m2, sum all member of sin_e*1360, all non-zero value will take part in calculation
    RadTot = RadTot*0.0036 ## Convert W.hr/m2 to MJ/m2 = Theoritical Hourly Radiation
    Tt = rad/RadTot  ## Compare theoritical Radiation (RadTot) to measured total radiation(rad)
    cl = 2.33 - 3.33*Tt    ## Fraction cloud cover  12.4
    if cl<0:
        cl = 0
    if cl>1:
        cl = 1
    sin_e2 = sin_d * sin_la + cos_d * cos_la * np.cos(0.2618*(hour-12))
    if sin_e2>0:
        Sw_in = 1360*Tt*sin_e2 ## 12.7, Sw = ShortWave(W/m2)
    else: 
        Sw_in = 0          
    e_sat, de_sat = SVP(Ta)    
    eps_a  = 1.24* np.power((e_sat*RH/Ta),(1/7.)) ## SisPAT
    eps_ac = (1-0.84*cl)*eps_a + 0.84*cl      ## 12.3, atmospheric emissivity with cloud, if cl = 0 : eps_ac = eps_a
    Lw_in =  eps_ac * sb * np.power((Ta+273),4)   ## 5.6e-08 = Stefan Boltzman
    return Sw_in, Lw_in

def soil_Rn(Sw_in, Lw_in, LAI, alb_s, eps_s, Ts, Ta):
    tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
    Sw_s = tau_c * (1-alb_s) * Sw_in
    ##Lw_s = tau_c * (Lw_in - eps_s * sb * (mth.pow((Ta+273), 4) - mth.pow((Ts+273), 3)*(Ts-Ta))) ## Linearize Ts
    Lw_s = tau_c * eps_s *(Lw_in - sb * np.power((Ts+273), 4)) ##canopy or soil surface always compared to air
    Rn_s = Sw_s + Lw_s
    return Rn_s

def canopy_Rn(Sw_in, Lw_in, LAI, alb_c, eps_c, Tc, Ta):
    tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
    Sw_c = (1-tau_c) * (1-alb_c) * Sw_in
    ##Lw_c = (1-tau_c) * (Lw_in - eps_c * sb * (mth.pow((Ta+273), 4) - mth.pow((Tc+273), 3)*(Tc-Ta))) ## Linearize Tc
    Lw_c = (1-tau_c) * eps_c *(Lw_in-  sb * np.power((Tc+273), 4)) ##canopy or soil surface always compared to air
    Rn_c = Sw_c + Lw_c
    return Rn_c
###########################################################################

time = np.arange(0., 24., 0.1666666) ##create time array from 0 to 24 (1 day)
Ta = np.zeros(time.size) ##allocating space for Ta array
Tc = np.zeros(time.size) ##allocating space for Ta array
Ts = np.zeros(time.size) ##allocating space for Ta array
RH = np.zeros(time.size) ##allocating space for RH array
Sw_in = np.zeros(time.size) ##allocating space for Sw_in array
Lw_in = np.zeros(time.size) ##allocating space for Lw_in array
Rn_s = np.zeros(time.size) ##allocating space for Rn_s array
Rn_c = np.zeros(time.size) ##allocating space for Rn_c array

##time looping, start all process here
for j in xrange(time.size):
    Ta[j] = fourier(deltaTa, j, zTa, aTa0, aTa, bTa) ##fill array of Ta in 24 hours
    Tc[j] = fourier(deltaTc, j, zTc, aTc0, aTc, bTc) ##fill array of Ta in 24 hours
    Ts[j] = fourier(deltaTs, j, zTs, aTs0, aTs, bTs) ##fill array of Ta in 24 hours
    RH[j] = fourier(deltaRH, j, zRH, aRH0, aRH, bRH)/100 ##fill array of RH in 24 hours
    Sw_in[j], Lw_in[j] = Radiation(daystrt, time[j], rad, Ta[j], RH[j])      

for j in xrange(time.size):
    ##Calculate radiation in canopy surface and soil surface
    Rn_s[j] = soil_Rn(Sw_in[j], Lw_in[j], LAI, alb_s, eps_s, Ts[j], Ta[j])
    Rn_c[j] = canopy_Rn(Sw_in[j], Lw_in[j], LAI, alb_c, eps_c, Tc[j], Ta[j])

## Write value to spreedsheet. an example 
j = 0
for cell in pyo('2-Layer.f6:f150'):
    cell.value = Rn_s[j]
    j = j + 1

plt.figure(1)
## Plot Global Radiation Sw_in
plt.subplot(211)
plt.plot(time, Sw_in, 'b-')
plt.title('Shortwave Incoming Radiation')
plt.xlabel('time (hours)')
plt.ylabel('Incoming Radiation (W/m^2)')
plt.subplot(212)
plt.plot(time, Lw_in, 'r-')
plt.title('Longwave Incoming Radiation')
plt.xlabel('time (hours)')
plt.ylabel('Incoming Radiation (W/m^2)')
plt.figure(2)
## Plot Soil Net Radiation
plt.subplot(211)
plt.plot(time, Rn_s, 'r--')
plt.title('Net Radiation in Soil')
plt.xlabel('time (hours)')
plt.ylabel('Soil Net Radiation (W/m^2)')
## Plot Canopy Net Radiation
plt.subplot(212)
plt.plot(time, Rn_c, 'r--')
plt.title('Net Radiation in Canopy')
plt.xlabel('time (hours)')
plt.ylabel('Canopy Net Radiation (W/m^2)')

plt.show()

</pre> 
