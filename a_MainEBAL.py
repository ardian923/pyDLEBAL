#! /usr/bin/env python
## Soil-Plant-Atmospheric Continuum calculation emphasizing on surface energy balance
## Developed initially by Ardiansyah (UNSOED), http://ardiansyah.net
# ask for permission from me : ardi.plj@gmail.com (cc : ard@unsoed.ac.id) if you want to use this. It's free

# import all module from file, reading all parameters in spreadsheet
import os as os ##operating system function
from oosheet import OOSheet as pyo ##open office sheet
import b_doubleLayerEB as et
import c_WatEvapFlow as wef
import d_HeatFlow as hf
import e_FuncLibrary as fl


## Simulation parameters
dt = pyo('Inputs.c139').value; dt = float(dt) ## timestep
inittime = pyo('Inputs.c140').value
endtime = pyo('Inputs.c141').value
## Define numerical variables    
m = n - 2
z = hf.linGrid(depth, n , dz1) ## linear grid
## plot distance from surface
for col in xrange(1, m+2, 1):
    pyo('Soil.aa6')[0][col].value = z[col]

    Tb =  ## Initial soil temperature for each layer 
    Tn = initSoilTemp(Ta, Tb) ## Soil temperature profile 
    AA = np.zeros(n); BB = np.zeros(n); CC = np.zeros(n); DD = np.zeros(n) ## tridiagonal matrix coefficient

def JacobianMatric
    ## Prepare Matrix
    global A
    A = np.zeros(25).reshape(5,5)

    ## Jacobian Marix Calculation

    # Derivative over p_s
    dRn_s_p_s = 0                       ; dRn_c_p_s = 0
    dH_s_p_s  = 0                       ; dH_c_p_s  = 0       ; dH_a_p_s  = 0
    dLE_s_p_s = s.dLEsdps(p_s, Ts, e_o) ; dLE_c_p_s = 0       ; dLE_a_p_s = 0
    dG_s_p_s  = 0
    A[0 , 0] = dF_1_p_s  = dRn_s_p_s - dH_s_p_s - dLE_s_p_s - dG_s_p_s
    A[1 , 0] = dF_2_p_s  = dRn_c_p_s - dH_c_p_s - dLE_c_p_s
    A[2 , 0] = dF_3_p_s  = dH_a_p_s - dH_s_p_s - dH_c_p_s
    A[3 , 0] = dF_4_p_s  = dLE_a_p_s - dLE_s_p_s - dLE_c_p_s
    A[4 , 0] = dF_5_p_s  =

    # Derivative over e_o
    dRn_s_e_o = 0                       ; dRn_c_e_o = 0
    dH_s_e_o  = 0                       ; dH_c_e_o  = 0                 ; dH_a_e_o  = 0
    dLE_s_e_o = s.dLEsde_o(p_s, Ts, e_o); dLE_c_e_o = dLEcde_o(Tc, e_o) ; dLE_a_e_o = a.dLEade_o(Ta, e_o)
    dG_s_e_o  = 0
    A[0 , 1] = dF_1_e_o  = dRn_s_e_o - dH_s_e_o - dLE_s_e_o - dG_s_e_o
    A[1 , 1] = dF_2_e_o  = dRn_c_e_o - dH_c_e_o - dLE_c_e_o
    A[2 , 1] = dF_3_e_o  = dH_a_e_o - dH_s_e_o - dH_c_e_o
    A[3 , 1] = dF_4_e_o  = dLE_a_e_o - dLE_s_e_o - dLE_c_e_o
    A[4 , 1] = dF_5_e_o  =

    # Derivative over Ts
    dRn_s_Ts = dRnsdTs(LAI, alb_s, eps_s, Ts, Ta)                       ; dRn_c_Ts = 0
    dH_s_Ts  = s.dHsdTs(Ts, T0)         ; dH_c_Ts  = 0                  ; dH_a_Ts  = 0
    dLE_s_Ts = s.dLEsdTs(p_s, Ts, e_o)  ; dLE_c_Ts = 0                  ; dLE_a_Ts = 0
    dG_s_Ts  =
    A[0 , 2] = dF_1_Ts  = dRn_s_Ts - dH_s_Ts - dLE_s_Ts - dG_s_Ts
    A[1 , 2] = dF_2_Ts  = dRn_c_Ts - dH_c_Ts - dLE_c_Ts
    A[2 , 2] = dF_3_Ts  = dH_a_Ts - dH_s_Ts - dH_c_Ts
    A[3 , 2] = dF_4_Ts  = dLE_a_Ts - dLE_s_Ts - dLE_c_Ts
    A[4 , 2] = dF_5_Ts  =

    # Derivative over Tc
    dRn_s_Tc = 0                        ; dRn_c_Tc = p.dRncdTc(LAI, alb_c, eps_c, Tc, Ta)
    dH_s_Tc  = 0                        ; dH_c_Tc  = p.dHcdTc(Tc, T0)   ; dH_a_Tc  = 0
    dLE_s_Tc = 0                        ; dLE_c_Tc = p.dLEcdTc(Tc, e_o) ; dLE_a_Tc = 0
    dG_s_Tc  = 0
    A[0 , 3] = dF_1_Tc  = dRn_s_Tc - dH_s_Tc - dLE_s_Tc - dG_s_Tc
    A[1 , 3] = dF_2_Tc  = dRn_c_Tc - dH_c_Tc - dLE_c_Tc
    A[2 , 3] = dF_3_Tc  = dH_a_Tc - dH_s_Tc - dH_c_Tc
    A[3 , 3] = dF_4_Tc  = dLE_a_Tc - dLE_s_Tc - dLE_c_Tc
    A[4 , 3] = dF_5_Tc  =

    # Derivative over T0
    dRn_s_T0 = 0                        ; dRn_c_T0 = 0
    dH_s_T0  = s.dHsdT0(Ts, T0)         ; dH_c_T0  = p.dHcdT0(Tc, T0)   ; dH_a_T0  = a.dHadT0(T0, Ta)
    dLE_s_T0 = 0                        ; dLE_c_T0 = 0                  ; dLE_a_T0 = 0
    dG_s_T0  = 0
    A[0 , 4] = dF_1_T0  = dRn_s_T0 - dH_s_T0 - dLE_s_T0 - dG_s_T0
    A[1 , 4] = dF_2_T0  = dRn_c_T0 - dH_c_T0 - dLE_c_T0
    A[2 , 4] = dF_3_T0  = dH_a_T0 - dH_s_T0 - dH_c_T0
    A[3 , 4] = dF_4_T0  = dLE_a_T0 - dLE_s_T0 - dLE_c_T0
    A[4 , 4] = dF_5_T0  =
    return A

B = np.array([p_s, e_o, Ts, Tc, T0])
dB = np.array([dp_s, de_o, dTs, dTc, dT0])

time = 0
i = 0
## time looping
while (time <= endtime):
    
    daystrt = et.daystrt
    Ta = ## based on time, read data
    RH = ## based on time, read data
    rad= ## based on time, read data
    alb_s = et.alb_s
    alb_c = et.alb_c
    eps_s = et.eps_s
    eps_c = et.eps_c
    LAI = et.LAI


    ## Energy balance calculation (with or without vegetation all depends to LAI)
    s = et.soil(); p = et.plant(); a = et.atmosphere()

    Sw_in, Lw_in = a.Radiation(daystrt, time, rad, Ta, RH) ## hour, rad, Ta, RH in array
    Lw = (Lw_in - eps_s * sb * np.power((Ta+273),4) ) ## net Long wave radiation for isothermal


    ##Set initial guess for unknown parameters
    e_o = ## based on time, data or initial guess    
    T0  = ## based on time, data or initial guess
    Tc  = ## based on time, data or initial guess
    Ts  = ## based on time, data or initial guess
    p_s = ## based on time, data or initial guess    

    ## run soil heat flow
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

    ##Calculate radiation in canopy surface and soil surface
    Rn_s = s.soil_Rn(Sw_in, Lw_in, LAI, alb_s, eps_s, Ts, Ta)
    Rn_c = p.canopy_Rn(Sw_in, Lw_in, LAI, alb_c, eps_c, Tc, Ta)
    H_s = s.soil_H(Ts, T0)
    H_c = p.canopy_H(Tc, T0) 
    H_a = a.H(T0, Ta) 
    LE_s = s.soil_LE(Ts, e_o)
    LE_c = p.canopy_LE(Tc, e_o)
    LE_a = a.LE(Ta, e_o)
    G_s = Ghf

    ##Call all two-layer model energy balance equations
    F_1 = Rn_s - H_s - LE_s - G_s  ##f(p_s, Ts, T0, e_o)
    F_2 = Rn_c - H_c  - LE_c       ##f(Tc, T0, e_o)
    F_3 = H_a - H_s - H_c          ##f(Ts, Tc, T0)
    F_4 = LE_a - LE_s - LE_c       ##f(p_s, Ts, Tc, e_o) 
    F_5 = (Flux-ETs)+B1+C2-D1;
    
    ## Solve for ps, Ts, T0, Tc, e_o
    A = JacobianMatric
    ## Update all node of T and p in soil 

    ## write result to spreedsheet


    # update time step
    time = time + dt/3600.
    i = i + 1

##################### ==NOTES== ################################
















