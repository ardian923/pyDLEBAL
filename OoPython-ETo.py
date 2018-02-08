#! /usr/bin/env python
## Soil-Plant-Atmospheric Continuum calculation emphasizing on surface energy balance
## Developed initially by Ardiansyah (UNSOED), http://ardiansyah.net

##ardiansyah@AL-FATIH-II:~/Desktop/oosheet-1.0$ oocalc OoPython-ETo.ods -accept="socket,host=localhost,port=2002;urp;StarOffice.ServiceManager"
##python

import os as os ##operating system function
from oosheet import OOSheet as pyo ##open office sheet
import numpy as np ##numerical python library
import matplotlib.pyplot as plt
import scipy as sc ##scientific python library
import math as mth ##matemathical operation library
 
##obtain all input from spreedsheet
##consult spreedsheet for parameters note, symbol, unit, etc
##Boundary Layer Parameters
LAI = pyo('Inputs.c8').value
LAImax = pyo('Inputs.c9').value 
u = pyo('Inputs.c10').value
ha = pyo('Inputs.c11').value
hc = pyo('Inputs.c12').value
alt = pyo('Inputs.c13').value

##Initial Resistances
rss = pyo('Inputs.c17').value
rsa = pyo('Inputs.c18').value
rcs = pyo('Inputs.c19').value
rca = pyo('Inputs.c20').value
raa = pyo('Inputs.c21').value

##Radiation Parameters
alb_s = pyo('Inputs.c25').value
alb_c = pyo('Inputs.c26').value
eps_s = pyo('Inputs.c27').value
eps_c = pyo('Inputs.c28').value
deg = pyo('Inputs.c30').value
mnt = pyo('Inputs.c31').value
sec = pyo('Inputs.c32').value
hem = pyo('Inputs.c33').string
rad = pyo('Inputs.c34').value
daystrt = pyo('Inputs.c35').value

##Initial Condition
e_o = pyo('Inputs.c38').value
T0 = pyo('Inputs.c39').value
Tc = pyo('Inputs.c40').value
Ts = pyo('Inputs.c41').value
p_s = pyo('Inputs.c42').value
dz1 = pyo('Inputs.c43').value

##Constants
gr = pyo('Inputs.c48').value
R = pyo('Inputs.c49').value
Rsp = pyo('Inputs.c50').value
Ma = pyo('Inputs.c51').value
rho_a = pyo('Inputs.c52').value  #function available
Cp = pyo('Inputs.c53').value
Mw = pyo('Inputs.c54').value
rho_w = pyo('Inputs.c55').value
eps = pyo('Inputs.c56').value
lamda = pyo('Inputs.c57').value  #function available
Dv = pyo('Inputs.c58').value
gamma = pyo('Inputs.c59').value  #function available
Vk = pyo('Inputs.c60').value
p_atm = pyo('Inputs.c61').value  #function available
sb = pyo('Inputs.c62').value

##Constant for fourier RH and Ta
delta = pyo('Inputs.c66').value
##Constant for fourier Ta
aT0 = pyo('Inputs.c68').value
aT = np.array([pyo('Inputs.c69').value, pyo('Inputs.c70').value, pyo('Inputs.c71').value, pyo('Inputs.c72').value])
bT = np.array([pyo('Inputs.c73').value, pyo('Inputs.c74').value, pyo('Inputs.c75').value, pyo('Inputs.c76').value])
##Constant for fourier RH
aRH0 = pyo('Inputs.c104').value
aRH = np.array([pyo('Inputs.c105').value, pyo('Inputs.c106').value, pyo('Inputs.c107').value, pyo('Inputs.c108').value])
bRH = np.array([pyo('Inputs.c109').value, pyo('Inputs.c110').value, pyo('Inputs.c111').value, pyo('Inputs.c112').value])

os.system('clear') #clear screen on linux

def AtmPressure(alt): ## p_atm, in kPa, altitude in meter
    AtmPressure = 101.3 * np.exp(-alt/8200)
    return AtmPressure

def AtmDensity(Pa, Ta): ## rho_a, Atmospheric density (density of air), kg/m3
    AtmDensity = (1000 * Pa)/(1.01 * Rsp * (Ta+273.16))   ## Pa = atmospheric pressure, Rsp = 287 J/kg/K
    return AtmDensity

def SVP(T): ## Saturated Vapor Pressure
        e_sat  = 0.611*np.exp(17.502*T/(T+240.97)) ## 3.8, Campbell 1998 (kPa)
        de_sat = 17.502*240.97*e_sat/np.sqrt(240.97+T)   ## 3.9, Campbell 1998
        return e_sat, de_sat ##after return we repeat the name of variable/s so that we can assign it later to another new variable

def AVP(T, Twb, alt): ## Actual Vapor Pressure
    e_w, de_w = SVP(Twb)
    Pa = AtmPressure(alt) * 10 ## in hPa
    e_act = e_w - Pa * (T - Twb) * 0.00066 * (1 + (0.00115 * Twb))
    return e_act

def LatentHeatVaporize(T): ##lambda, Latent heat of vaporization J/kg 
    LatentHeatVaporize = (45144 - 48*T)/(Mw)
    return LatentHeatVaporize

def VaporDiffus(patm,T): ##Vapor Difussivity (Dv) in air
    VaporDiffus = 0.0000214*(101.3/patm)*np.power((T/273.16),1.75)       ##7.7, Campbell 1998
    return VaporDiffus

def PsycConstant(Pa, T): ## gamma, Psychrometric Constant (kPa/K)
    lamda = LatentHeatVaporize(T)
    PsycConstant = (Cp * Pa) / (eps * lamda)   ## Pa = atmospheric pressure
    return PsycConstant

def fourier(delta, j, z, a0, a, b): ##create fourier function for air temperature and air humidity
    i = np.arange(1, 5) ##create i array for sigma 
    w = 2 * np.pi / 24
    sigma = np.exp(-z*(i*w/(2*delta))**0.5) * (a * np.cos(i*w*time[j] - z*(i*w/(2*delta))**0.5)) + (b * np.sin(i*w*time[j] - z*(i*w/(2*delta))**0.5))
    fourier = a0 + sum(sigma)
    return fourier

def ConvertLattitude(deg, mnt, sec, hem):
    LattRadian = deg+(mnt/60.0)+(sec/3600.0)
    LattRadian = (np.pi/180)*LattRadian
    if hem == 'south':
        LattRadian =-LattRadian 
    else:
        LattRadian = LattRadian
    return LattRadian

## for Penman-Monteith ET calculation
def Simple_Go(time, Rn):
    ## calculate Go as fraction of Rn
    i = 0
    for x in time:
        if ((x-6)>=0 and (x-18)<0):
            Go[i] = 0.1 * Rn[i]   ## daytime Go for tall crop (>0.5 m) crop = 0.5, for short 0.1*Rn
        else:
            Go[i] = 0.04 * Rn[i]  ## nightime Go for tall crop (>0.5 m) crop = 0.2, for short 0.04*Rn
        i = i + 1
    return Go

def PenmanMonteith_ETo(u, Ta, Twb, ha, Rn, Go):
    lamda = LatentHeatVaporize(Ta)
    e_sat, de_sat = SVP(Ta)   
    e_act = AVP(Ta, Twb, alt)
    VPD = e_sat - e_act   
    Pa = AtmPressure(alt)
    gamma = PsycConstant(Pa, Ta)
    rho_a = AtmDensity(Pa, Ta)     
    u2 = u * (4.87)/(np.log(67.8*ha - 5.42))
    raa = 208./u2 				   ## aerodynamic resistance (s/m)
    rca = 30 ## daytime = 50 (short plant <0.5 m), 30 (tall plant > 0.5), nighttime = 200 (s/m)
    PenmanMonteith_ETo = (1/lamda) * (de_sat * (Rn - Go) +  (VPD * rho_a * Cp)/raa)/(de_sat + gamma * (1+(rca/raa)))
    return PenmanMonteith_ETo
## end Penman-Monteith ET calculation

class soil:
    'physical process in soil surface and below'
    def rss(self, dz1, por): ##soil resistance for evaporation, por = porosity
        rss = dz1/(0.66 * Dv * por) ##Dv = vapor diffusivity (m2/s), dz1 wet 0.0001
        ##give porosity = 0.5 if it is not calculatedin program
        return rss ##rss can't be zero if we want to put it as denominator

    def rsa(self, hc):   ##from Teh, eq. 4.68
        na = 2        ##2 - 3 (Monteith, 1975), attenuation coefficient for eddy diffusivity
        z_o   = 0.13*hc
        z_so  = 0.004   ##Table 4.1, Teh, for bare soil
        d  = 0.778*hc 
        ustar = Vk*u/(np.log((ha-d)/z_o))
        Kcan = Vk*ustar*hc 
        u1   = np.exp(-na*(z_so/hc)) 
        u2   = np.exp(-na*((z_o+d)/hc)) 
        rsa = ((hc*np.exp(na))/(na*Kcan))*(u1-u2)
        return rsa

    def soil_RH(self, p_soil, T_soil):
        ##p = (1/10.2)*p; //convert from cmH2O to kPa
        soil_RH = np.exp(Mw * p_soil / (R * (T_soil+273)))        ##eq 9.3 Campbell 1985
        return soil_RH

    def soil_LE(self, p_s, Ts, e_o): ## self. is calling function within class
        e_s, de_s = SVP(Ts)
        soil_LE = (rho_a * Cp/gamma)*(e_s* self.soil_RH(p_s, Ts) - e_o) /(self.rss(dz1, por) + self.rsa(hc)) 
        return soil_LE

    def soil_H(self, Ts, T0):
        soil_H = (rho_a * Cp) * (Ts - T0)/self.rsa(hc)
        return soil_H

    def soil_Rn(self, Sw_in, Lw_in, LAI, alb_s, eps_s, Ts, Ta):
        tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
        Sw_s = tau_c * (1-alb_s) * Sw_in
        ##Lw_s = tau_c * (Lw_in - eps_s * sb * (mth.pow((Ta+273), 4) - mth.pow((Ts+273), 3)*(Ts-Ta))) ## Linearize Ts
        Lw_s = tau_c * eps_s *(Lw_in - sb * np.power((Ts+273), 4)) ##canopy or soil surface always compared to air
        Rn_s = Sw_s + Lw_s
        return Rn_s
    ## derivative some function
    def dLEsdps(self, p_s, Ts, e_o): ## self. is calling function within class
        e_s, de_s = SVP(Ts)
        dLEsdRH = (rho_a * Cp/gamma)*(e_s) /(self.rss(dz1, por) + self.rsa(hc)) 
        dLEsdps = dLEsdRH * self.soil_RH(p_s, Ts) * (Mw / (R * (Ts+273)))     
        return dLEsdps

    def dLEsdTs(self, p_s, Ts, e_o): ## self. is calling function within class
        e_s, de_s = SVP(Ts)
        dRHdTs  = self.soil_RH(p_s, Ts) * (- Mw * p_soil / R )  * np.power((Ts+273),2)
        dLEsdTs = (rho_a * Cp/gamma)*(de_s*self.soil_RH(p_s, Ts) + e_s*dRHdTs) /(self.rss(dz1, por) + self.rsa(hc)) 
        return dLEsdTs

    def dLEsde_o(self, p_s, Ts, e_o): ## self. is calling function within class
        e_s, de_s = SVP(Ts)
        dLEsde_o = (rho_a * Cp/gamma)*(- e_o) /(self.rss(dz1, por) + self.rsa(hc)) 
        return dLEsde_o

    def dHsdTs(self, Ts, T0):
        dHsdTs = (rho_a * Cp) /self.rsa(hc)
        return dHsdTs

    def dHsdT0(self, Ts, T0):
        dHsdT0 = -(rho_a * Cp)/self.rsa(hc)
        return dHsdT0

    def dRnsdTs(self, LAI, alb_s, eps_s, Ts, Ta):
        tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
        dRnsdTs = tau_c * eps_s *(- sb * 4 * np.power((Ts+273), 3))
        return dRnsdTs

class plant:
    'physical process in soil plant canopy'
#    def rcs(self, LAI):  ##canopy surface resistance
#        ##when leaf water potential decreases below critical,stomatal resistance increase rapidly
#        sp = 10
#        rsto0 = 100 ##rsto0 = minimum resistance of stomata, s/m
#        rsto = rsto0*(1+np.power((pl,pc),sp)) ##sp = stomata critical constant = 10 (arbitary)
#        rcs = rsto/LAI
#        return rcs ################################################################### ATTENTION!! THIS STILL FRAGILE but simple #################
    def rcs (self, rsto, LAI, LAImax): ## Teh 4.80
        if LAI <= LAImax:
            rcs = rsto/LAI  
        else:
            rcs = rsto/(0.5*LAImax)
        return rcs

    def rca(self, LAI): ##canopy resistane or bulk stomatal resistance
        rca = 88.816*np.exp(0.4192*LAI) ##Campbell-Soil physics with Basic
        return rca

    def canopy_LE(self,Tc, e_o):
        e_sc, de_sc = SVP(Tc)
        canopy_LE = (rho_a * Cp/gamma) * (e_sc - e_o) /(self.rcs(LAI) + self.rca(LAI))
        return canopy_LE

    def canopy_H(self,Tc, T0):
        canopy_H = (rho_a * Cp) * (Tc - T0)/self.rca(LAI)
        return canopy_H

    def canopy_Rn(self, Sw_in, Lw_in, LAI, alb_c, eps_c, Tc, Ta):
        tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
        Sw_c = (1-tau_c) * (1-alb_c) * Sw_in
        ##Lw_c = (1-tau_c) * (Lw_in - eps_c * sb * (mth.pow((Ta+273), 4) - mth.pow((Tc+273), 3)*(Tc-Ta))) ## Linearize Tc
        Lw_c = (1-tau_c) * eps_c *(Lw_in-  sb * np.power((Tc+273), 4)) ##canopy or soil surface always compared to air
        Rn_c = Sw_c + Lw_c
        return Rn_c
    ## derivative some function
    def dLEcdTc(self,Tc, e_o):
        e_sc, de_sc = SVP(Tc)
        dLEcdTc = de_sc * (rho_a * Cp/gamma) /(self.rcs(LAI) + self.rca(LAI))
        return dLEcdTc

    def dLEcde_o(self,Tc, e_o):
        e_sc, de_sc = SVP(Tc)
        dLEcde_o = - (rho_a * Cp/gamma)/(self.rcs(LAI) + self.rca(LAI))
        return dLEcde_o

    def dHcdTc(self,Tc, T0):
        dHcdTc = (rho_a * Cp) /self.rca(LAI)
        return dHcdTc

    def dHcdT0(self,Tc, T0):
        dHcdT0 = - (rho_a * Cp) /self.rca(LAI)
        return dHcdT0

    def dRncdTc(self, LAI, alb_c, eps_c, Tc, Ta):
        tau_c = np.exp(-0.4*LAI)   ## tau_c = plant transmittance coefficient, no vegetation, tau_c = 1     
        dRncdTc = (1-tau_c) * eps_c * (- sb * 4 * np.power((Tc+273), 3))
        return dRncdTc


class atmospher:
    'physical process above canopy'
    def raa(self,ha, u, hc, Tsurf, Ta):  ## Boundary layer resistance, Tsurf = T0(vegetated) or Ts (bare soil)
        pm = 0
        ph = 0
        d  = 0.778*hc
        zm = 0.13*hc
        zh = 0.2*zm
        a = np.arange(1, 4, 1)
        for i in xrange(a.size):
            ustar = u * Vk/(np.log((ha-d+zm)/zm)+pm)## u, mean wind speed  m/s
            raa = (np.log((ha-d+zh)/zh)+ph)/(Vk*ustar)
            K_h0 = (rho_a*Cp)/raa  ##K_h0: J/s.m2.K or W/m2.K
            sp = -Vk*ha*gr*K_h0*(Tsurf-Ta)/((rho_a*Cp)*(Ta+273)*np.power(ustar,3)) ##Soil physic with basic- 12.14, stability parameter
            if sp > 0:
                ph = 4.7*sp
                pm = ph
            else:
                ph = -2*np.log((1+np.sqrt(1-16*sp))/2)
                pm = 0.6*ph
        return raa

    def LE(self, Ta, e_o):
        e_sa, de_sa = SVP(Ta)
        e_a = e_sa * RH
        LE = (rho_a * Cp/gamma) * (e_o - e_a) /(self.rcs(LAI) + self.rca(LAI))
        return LE

    def H(self, T0, Ta):
        H = (rho_a * Cp) * (T0 - Ta)/self.raa(ha, u, hc, Tsurf, Ta)
        return H

    def Radiation(self, daystrt, hour, rad, Ta, RH): ## Calculate Radiation using input hour
        ## rad is total daily measured global radiation (MJ/m2)
        la = ConvertLattitude(deg, mnt, sec, hem)
        sin_la = np.sin(la) ## sin lattitude
        cos_la = np.cos(la) ## cos lattitude
        sin_d = 0.3985*np.sin(4.869+0.0172*daystrt+0.03345*np.sin(6.224+0.0172*daystrt))  ## sin declination 12.6, Campbell
        cos_d = np.sqrt(1-np.power(sin_d, 2))
        ## Prepare to calculate Total Radiation
        hrs = np.arange(0., 24., 0.1666666) ## calculate for every 10 minutes (600 sec) was :  np.arange(1., 25., 1)
        sin_e = np.zeros(hrs.size)
        sin_e = sin_d*sin_la + cos_d*cos_la*np.cos(0.2618*(hrs-12))  ## Solar elevation angle 12.8, returns array of sin_e
        RadTot = 0
#        for i in xrange(hrs.size):
#            if sin_e[i]>0:
#                RadTot = RadTot + 1360*sin_e[i] ## W-hr/m2, sum all member of sin_e*1360, all non-zero value will take part in calculation
#        RadTot = RadTot*(600./1e6) ## was : RadTot = RadTot*0.0036 ## Convert W.hr/m2 to MJ/m2 = Theoritical Hourly Radiation
        RadTot = 1360*sin_e
        Tt = rad/RadTot  ## Compare theoritical Radiation (RadTot) to measured total radiation(rad), if both are arrays, Tt would be an array
        cl = 2.33 - 3.33*Tt    ## Fraction cloud cover  12.4
#        if cl<0:
#            cl = 0
#        if cl>1:
#            cl = 1
        cl[cl < 0] = 0 ## replace negative value with zero
        cl[cl > 1] = 1 ## replace value larger than 1 to 1
#        sin_e2 = sin_d * sin_la + cos_d * cos_la * np.cos(0.2618*(hour-12))
#        if sin_e2>0:
#            Sw_in = 1360*Tt*sin_e2 ## 12.7, Sw = ShortWave(W/m2)
#        else: 
#            Sw_in = 0          
        sin_e2 = np.zeros(hour.size) ## need to allocate array space for sin_e2 since it doesn't have array variable/s
        sin_e2 = sin_d * sin_la + cos_d * cos_la * np.cos(0.2618*(hour-12))
        sin_e2[sin_e2 <= 0] = 0
        Sw_in = 1360*Tt*sin_e2
        e_sat, de_sat = SVP(Ta)    
        eps_a  = 1.24* np.power((e_sat*RH/Ta),(1/7.)) ## SisPAT
        eps_ac = (1-0.84*cl)*eps_a + 0.84*cl      ## 12.3, atmospheric emissivity with cloud, if cl = 0 : eps_ac = eps_a
        Lw_in =  eps_ac * sb * np.power((Ta+273),4)   ## 5.6e-08 = Stefan Boltzman
        return Sw_in, Lw_in

    ## derivative some function
    def dLEade_o(self, Ta, e_o):
        dLEade_o = (rho_a * Cp/gamma) * /(self.rcs(LAI) + self.rca(LAI))
        return dLEade_o

    def dHadT0(self, T0, Ta):
        dHadT0 = (rho_a * Cp) / self.raa(ha, u, hc, Tsurf, Ta)
        return dHadT0




###############################################################
if __name__ == '__main__':  ## Run simulation as standalone program
    ## Main Program
    ##assign s-p-a to each class
    s = soil()
    p = plant()
    a = atmospher()

    ##allocating space for array variable
    time = np.arange(0., 24., 0.1666666) ##create time array from 0 to 24 (1 day)
    ## sensor measurement input or Fourier input or data input
    Ta = np.zeros(time.size) ##allocating space for Ta array
    Twb = np.zeros(time.size) ##allocating space for Ta wet bulb array
    rad = np.zeros(time.size) ##allocating space for incoming radiation
    RH = np.zeros(time.size) ##allocating space for RH array, calculated from Ta wet bulb
    ## estimated or calculated during runtime
    Tc = np.zeros(time.size) ##allocating space for Tc array
    Ts = np.zeros(time.size) ##allocating space for Ts array
    Sw_in = np.zeros(time.size) ##allocating space for Sw_in array
    Lw_in = np.zeros(time.size) ##allocating space for Lw_in array
    Rn_s = np.zeros(time.size) ##allocating space for Rn_s array
    Rn_c = np.zeros(time.size) ##allocating space for Rn_c array
    Go = np.zeros(time.size) ##allocating space for Go array
    z = 0 ## dumping depth?

    # obtain data value 
    i = 0
    for row in pyo('1-Layer.h13:h157'):
        Ta[i]  = pyo('1-Layer.h13:h157')[i]['H'].value
        Twb[i] = pyo('1-Layer.h13:h157')[i]['I'].value
        rad[i] = pyo('1-Layer.h13:h157')[i]['J'].value 
        i = i + 1

    ## Penman-Monteith procedure, try without looping
    e_sat, de_sat = SVP(Ta)   
    e_act = AVP(Ta, Twb, alt)
    RH = e_act/e_sat  

    Sw_in, Lw_in = a.Radiation(daystrt, time, rad, Ta, RH) ## hour, rad, Ta, RH in array
    Lw = (Lw_in - eps_s * sb * np.power((Ta+273),4) ) ## net Long wave radiation for isothermal
    Rn = 0.77*Sw_in + Lw ## Rn is array of isothermal net radiation
    Go = Simple_Go(time, Rn)
    ETo = PenmanMonteith_ETo(u, Ta, Twb, ha, Rn, Go)

##time looping, start all process here, looping ahould be avoided in python, but this optimization procedure needs it
#for j in xrange(time.size):
#    ##Ta[j] = fourier(delta, j, z, aT0, aT, bT) ##fill array of Ta in 24 hours
#    ##RH[j] = fourier(delta, j, z, aRH0, aRH, bRH) ##fill array of RH in 24 hours
#    Ta[j] = {read from spreadsheet}
#    RH[j] = {read from spreadsheet}
#    Sw_in[j], Lw_in[j] = a.Radiation(daystrt, time[j], rad, Ta[j], RH[j])  
#
#    ##Set initial guess for unknown parameters    
#    
#    ##Calculate radiation in canopy surface and soil surface
#    Rn_s[j] = s.soil_Rn(Sw_in[j], Lw_in[j], LAI, alb_s, eps_s, Ts[j], Ta[j])
#    Rn_c[j] = p.canopy_Rn(Sw_in[j], Lw_in[j], LAI, alb_c, eps_c, Tc[j], Ta[j])
#    H_s[j] = s.soil_H(Ts[j], T0[j])
#    H_c[j] = p.canopy_H(Tc[j], T0[j]) 
#    H_a[j] = a.H(T0[j], Ta[j]) 
#    LE_s[j] = s.soil_LE(Ts[j], e_o[j])
#    LE_c[j] = p.canopy_LE(Tc[j], e_o[j])
#    LE_a[j] = a.LE(Ta[j], e_o[j])
#   
#    
#    ##Call all two-layer model energy balance equations
#    F_1 := Rns - s.soil_H(Ts, T0) - s.soil_LE(Ts, e_o) - s.soil_G  ##f(p_s, Ts, T0, e_o)
#    F_2 := Rnc - p.canopy_H(Tc, T0) - p.canopy_LE(Tc, e_o)         ##f(Tc, T0, e_o)
#    F_3 := a.H - s.soil_H(Ts, T0) - p.canopy_H(Tc, T0)             ##f(Ts, T0)
#    F_4 := a.LE - s.soil_LE - p.canopy_LE(Tc, e_o)                 ##f(p_s, Ts, Tc, e_o) 
#    F_5 := (Flux-ETs)+B1+C2-D1;
#    F_6 := (Rns + Rnc) - s.soil_G - (s.soil_H + p.canopy_H) - (s.soil_LE + p.canopy_LE)
### Solve for ps, Ts, T0, Tc, e_o
###Update all node of T and p in soil 

    ##write result to spreedsheet
    i = 0
    for row in pyo('1-Layer.h13:h157'):
        pyo('1-Layer.h13:h157')[i]['F'].value = ETo[i]
        i = i + 1

    # plot results
    plt.figure(1)
    plt.subplot(111)
    plt.plot(time, Rn, 'r-')
    #plt.plot(time, rad, 'b-')
    #plt.plot(time, Lw, 'g-')
    #plt.plot(time, ETo, 'g-')
    #plt.legend(('ETo','Air Temperature', 'Relative Humidity'), 'upper right', shadow=False)
    #plt.title('Penman-Monteith Evapotranspiration')
    plt.xlabel('Time (24 hr decimal)')
    #plt.ylabel('ETo (W/m2)')
    plt.show()



