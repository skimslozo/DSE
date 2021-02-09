#Hydrogen Sizing

def getEnergyMass(Pavg, Pto, Tmission, Tvtol):
    t1 = Tvtol       #discharge time [s]
    t2 = Tmission/3.  #recharge time [s]
    ttot = Tmission   #total mission time 
    T = t1/t2     #Ratio               
    P1 = Pto #Peak power [W]
    P2 = Pavg  #Base power [W]
    Pf = (P1*t1+P2*t2)/(t1+t2)
    F = P1/P2    #Ratio
    n = 3. #number of cycles
    ec = P1*3*20.       #energy consumed
    Espec = 350.    #specific energy [Wh/kg]
    HHV = 39400.*3600.
    rhomh2 = 0.142    #table3
    rhovh2 = 70.8    #table3
    
    
    Pavg = (P1*t1+P2*t2)/(t1+t2)  #average power
    Es = n*(P1-Pavg)*(t1/3600) #Energy discharged from battery
    Efc = Pf*((t1+t2)/3600)*n   #Energy discharged from fuel cell
    
    if Pavg > P2 + 0.5*(P1-P2):
        Limit = "charge"
    if Pavg < P2 + 0.5*(P1-P2):
        Limit = "discharge"    
    
    npec = Es*(1/0.9)*(1/0.8)   #Name plate energy capacity
    Ws = (npec/Espec)*1.15            #Weight of battery [Wh/kg]
    Vs = (npec/Espec)*0.001           #Volume of battery [Wh/L]
    
    
    E = 33300 #energy density Wh/kg
    Eff = 0.6 #Total fuel cell efficiency
    
    Wh = (((Pavg*(ttot/3600))/Eff)/E)*1.15
    Vh = Wh/rhovh2 #volume hydrogen storage system
    
    Wt = Wh/0.7+Wh
    rhoPEM = 1500. #W/kg
    vrhopem = 1200. #W/L
    
    Wpem = Pavg/rhoPEM
    Vpem = (Pavg/vrhopem)*0.001
    
    Wtot = (Wpem + Wt + Ws)*1.2
    Vtot = (Vpem + Vh + Vs + (Vpem + Vh + Vs)*0.2)
    
    return Wtot

#print("Mbattery = " + str(Ws),
#      "Mhydrogen = " + str(Wh),
#      "Mtank = " + str(Wt),
#      "Mpem = " + str(Wpem),
#      "Mtot = " + str(Wtot))
#
#print("Vbattery = " + str(Vs),
#      "Vhydrogen = " + str(Vh),
#      "Vpem = " + str(Vpem),
#      "Vtot = " + str(Vtot))
#




