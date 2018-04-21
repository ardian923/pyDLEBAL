# pyDoubleLayerEB : Double Layer Model for Energy Balance in Bare and Vegetated Surface

Water Flow in Soil due to Infiltration or Evaporation
-----------------------------------------------------
Water flow in soil when there is gradient in soil matric pressure. The gradient could occur through infiltration (water in soil surface) or evaporation flux (vapor movement in soil surface)

**Liquid Flux :**

![](./pyDoubleLayerEB_files/equation.png)
		

where :
![](./pyDoubleLayerEB_files/equation010.png); ![](./pyDoubleLayerEB_files/equation011.png) = 40.3 ; ![](./pyDoubleLayerEB_files/equation012.png) = 0.672 ; ![](./pyDoubleLayerEB_files/equation013.png) = 3.31. The unit of  ![](./pyDoubleLayerEB_files/equation020.png) is depends on the unit of ![](./pyDoubleLayerEB_files/equation021.png)
Water retention curve is given by equation :
![](./pyDoubleLayerEB_files/equation006.png), where ![](./pyDoubleLayerEB_files/equation007.png) is matric potential in ![](./pyDoubleLayerEB_files/equation008.png), ![](./pyDoubleLayerEB_files/equation014.png) = 0.343; ![](./pyDoubleLayerEB_files/equation015.png) = 0.04133; ![](./pyDoubleLayerEB_files/equation016.png) = 5.603; ![](./pyDoubleLayerEB_files/equation017.png) = 0.8215; ![](./pyDoubleLayerEB_files/equation018.png) = 0.08624; ![](./pyDoubleLayerEB_files/equation019.png) = 14.37
water capacity is calculated as  ![](./pyDoubleLayerEB_files/equation009.png)


**Vapor Flux :**

![](./pyDoubleLayerEB_files/equation001.png); (![](./pyDoubleLayerEB_files/equation033.png)) - - -> if divided by ![](./pyDoubleLayerEB_files/equation034.png) (![](./pyDoubleLayerEB_files/equation035.png)), resulting flux in ![](./pyDoubleLayerEB_files/equation036.png)



where :
![](./pyDoubleLayerEB_files/equation002.png) ; (![](./pyDoubleLayerEB_files/equation032.png))

![](./pyDoubleLayerEB_files/equation003.png); (![](./pyDoubleLayerEB_files/equation038.png))

![](./pyDoubleLayerEB_files/equation005.png) is soil relative humidity

![](./pyDoubleLayerEB_files/equation004.png); ![](./pyDoubleLayerEB_files/equation022.png) (![](./pyDoubleLayerEB_files/equation031.png) - - -> if divided by 1000, resulting ![](./pyDoubleLayerEB_files/equation039.png)); ![](./pyDoubleLayerEB_files/equation023.png)

![](./pyDoubleLayerEB_files/equation029.png) = vapor difussivity (0.000024 ![](./pyDoubleLayerEB_files/equation030.png)); ![](./pyDoubleLayerEB_files/equation025.png) = gas constant (8.3143 ![](./pyDoubleLayerEB_files/equation024.png)); ![](./pyDoubleLayerEB_files/equation026.png) = mass of a mole water (0.018 ![](./pyDoubleLayerEB_files/equation027.png)); ![](./pyDoubleLayerEB_files/equation028.png) = temperature (Kelvin); ![](./pyDoubleLayerEB_files/equation037.png) = air filled porosity

```python
    kvi  = 0.66 * phi_i * Dv * (Mw / (R * (Ti[i] + 273))) * rhov_i  ## Dv, Mw inputted in spreedsheet
    kvTi = 0.66 * phi_i * Dv * drhov_sati * eta * hbar_i
    kvn  = 0.66 * phi_n * Dv * (Mw / (R * (Tn[i] + 273))) * rhov_n
    kvTn = 0.66 * phi_n * Dv * drhov_satn * eta * hbar_n
```
	

Heat Flow in Soil
-----------------
Heat flow in soil occures due to the existance of heat source. Soil surface explosed to solar radiation is the heat source.

Single Layer Energy Balance
---------------------------
Single layer energy balance assumes soil and plant are 1-layer that release vapor to air. Both evaporation (E) and transpiration (T) calculated as ET. Penman-Monteith equation is used to draw the relation.

```python
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
```

Penman-Monteith equation take soil and canopy as isothermal. First term of the eqiauation considering soil surface that requires Net Radiatin (R_n) and Ground Heat Flux (G_o). Second term of equation is about vapor transport from canopy that depends on canppy resistance (rca) and aerodynamic resistance (raa)

Double Layer Energy Balance
---------------------------

