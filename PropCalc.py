#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 09:17:52 2019
@author: MatsKlijn
"""
from math import *
from matplotlib import *
import numpy as np 
from readPars import load_input
import pandas as pd

missiondat=load_input(sheet='Mission')
acdat=load_input(sheet='Aircraft parameters')
energdat=load_input(sheet='Energy storage')
miscdat=load_input(sheet='Misc')
noisedat=load_input(sheet='Noise')

MTOW        = acdat.MTOW                                                         #MTOW from parameter excel in kg


def propellerdesign(MTOW, CD0, test=None):
    
    FOM         = acdat.M                                                            #Figure of Merit from the helicopter lady
    rho         = missiondat.rhoSL                                                   #At sea level
    rhotrans    = missiondat.rhoTrans                                                #Density at transitioning altitude
    ROC_vtol    = missiondat.ROCvtol
    TW_to       = 1.1
    Ti          = acdat.Ti
    M_tmax_c      = noisedat.M_tmax_c                                                               #Set a maximum tip speed to 0.8Mach
    sos_c         = missiondat.Ccruise                                                              #Speed of sound at sea level in m/s                                                 #range of number of possible propellers
    r_c          =  missiondat.Hcruise                                                                #Distance from sound source to sound receiver
    MTOW_c        = MTOW   
    Pset=acdat.Pset
    mup=acdat.mup
    CD0=CD0
    rhoCruise=missiondat.rhoCruise
    Vcruise=missiondat.Vcruise
    CLmax=acdat.CLmax
    Vstall=0.5*Vcruise
    Mpl=missiondat.Mpl
    ecc=acdat.e
    A=acdat.A    
    PNL_A = 88.5+4
    PNL = PNL_A + 14                                                                 #dB, requirement for aircraft at close range
    noyA = 10**((PNL_A-26)/(33+1/3))
    noy  = 10**((PNL-26)/(33+1/3))                                                #maximum value of noy to achieve the requirements for PNL A-weighted
    #print(noyA)
    #print(noy)
    noy_range = np.arange(1,noy,0.07)
    N_motors500 = 1
    # =============================================================================
    # acceptable_noys = [[50, 114],[63, 113],[80, 111],[100, 109],[125, 108],
    #                    [160, 107],[200, 105], [250, 104],[315,103],[400,102],
    #                    [500, 102],[630, 102],[800,102], [1000,102], [1250,100],
    #                    [1600,96],[2000,94],[2500,92],[3150,91],[4000,91],
    #                    [5000,92],[6300,93],[8000,96],[10000,99]]                     #[freq, SPL]
    # =============================================================================
    
    acceptable_noys = [[50, 132],[63, 131],[80, 129],[100, 127],[125, 126],
                       [160, 125],[200, 123], [250, 122],[315,121],[400,120],
                       [500, 120],[630, 120],[800,120], [1000,120], [1250,118],
                       [1600,114],[2000,112],[2500,110],[3150,109],[4000,109],
                       [5000,110],[6300,111],[8000,114],[10000,117]]                     #[freq, SPL]
    
    N_blades    = np.arange(2,7,1)                                                   #constraint on number of blades. Min = 2, Max = 5
    parameter   = []                                                                 #To append in for loop next lines, [RPM,N_blades]
    M_tmax      = 0.8                                                                #Set a maximum tip speed to 0.8Mach
    sos         = 343.                                                               #Speed of sound at sea level in m/s
    N_prop      = np.arange(4,5,1)                                                   #range of number of possible propellers
    r           = 50.                                                                #Distance from sound source to sound receiver
    
    
    N_prop2 = []
    N_prop3 = []
    N_prop4 = []
    N_prop5 = []
    N_prop6 = []
    N_prop7 = []
    N_prop8 = []
    
    for i in range(len(N_blades)):
        #print(i)
        for l in range(len(acceptable_noys)):
            RPM     = (acceptable_noys[l][0]*60)/N_blades[i]                         #RPMs
            #print(N_blades[i])
            D_prop  = (M_tmax*60*sos)/(pi*RPM)                                       #b = diameter propeller
            #print(D_prop)
            if RPM > 1000 and RPM < 3000: 
                if D_prop < 5.5 and D_prop > 0.5:                                        #Define range of propeller diameter between 0.5 and 7.5m
                    for m in range(len(N_prop)):
                        P_mot               = e**(1/15.3*(acceptable_noys[l][1]-83.4+20*log(D_prop)-38.5*M_tmax+3*(N_blades[i]-2)-10*log(N_prop[m])+20*log(r)))    #Formula for SPL rewritten to get P_motor       
                        P_mot_per_propeller = P_mot/N_prop[m]
                        #if P_mot_per_propeller < 500:                                    #Maximum technology motor now is 500kW so only take shit below this value
                        DL          = MTOW*9.81/(pi*D_prop**2/4)#/N_prop[m]                    #Disk loading formula in N/m^2 per propeller
                        P_req_hover = TW_to*sqrt(DL/(2*rhotrans*Ti**3))/FOM*MTOW*9.81/1000    #Power required for hover
                        Vh          = sqrt(DL/(2*rhotrans))                          #Vh is induced velocity (or mass flow)
                        P_to        = P_req_hover*((0.5*ROC_vtol/Vh)+sqrt(0.25*((ROC_vtol/Vh)**2+1)))
                        if P_to/N_prop[m] <  280:   
                        #print(P_to)
                            if DL > 1000 and DL < 10000:                                   #Assumed disk loading between 100 and 1000 N/m^2 due to lower disk loading is more efficient for hover
                                if P_to < P_mot: 
                                    if N_prop[m] > 1 and N_prop[m] < 3:
                                        N_prop2.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 2 and N_prop[m] < 4:
                                        N_prop3.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 3 and N_prop[m] < 5:
                                        N_prop4.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 4 and N_prop[m] < 6:
                                        N_prop5.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 5 and N_prop[m] < 7:
                                        N_prop6.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 6 and N_prop[m] < 8:
                                        N_prop7.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)])
                                    if N_prop[m] > 7 and N_prop[m] < 9:
                                        N_prop8.append([N_prop[m],RPM,N_blades[i],acceptable_noys[l][1],round(D_prop,2),round(P_mot,2),round(DL,2),round(P_to,2)]) 
                                    parameter.append(['RPM',RPM,'#Blades',N_blades[i],'#Propellers',N_prop[m],'SPl',acceptable_noys[l][1],'Diameter Prop',round(D_prop,2),'P_mot',round(P_mot,2),'DL',round(DL,2),'P_to',round(P_to,2)])     #Appended in 'parameter' is [RPM,N_blades, SPL value, N_diameter, P_mot, Disk Loading, P_req, N_propellers] related to eachother
    #print(len(parameter))                                               
    # =============================================================================
    # RPMs            = []
    # blade_amount    = [] 
    # Diameters       = [] 
    # SPLs            = []
    # P_motors        = []
    # DLs             = []
    # P_takeoffs      = []
    # 
    # for i in range(len(N_prop2)):
    #     RPMs.append(N_prop2[i][1])
    #     blade_amount.append(N_prop2[i][2])
    #     SPLs.append(N_prop2[i][3])
    #     Diameters.append(N_prop2[i][4])
    #     P_motors.append(N_prop2[i][5])
    #     DLs.append(N_prop2[i][6])
    #     P_takeoffs.append(N_prop2[i][7])
    # 
    # print('2 propellers \n'
    #       'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #       'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #       'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #       'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #       'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #       'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #       'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')                
    # =============================================================================
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #           
    #for i in range(len(N_prop3)):
    #    RPMs.append(N_prop3[i][1])
    #    blade_amount.append(N_prop3[i][2])
    #    SPLs.append(N_prop3[i][3])
    #    Diameters.append(N_prop3[i][4])
    #    P_motors.append(N_prop3[i][5])
    #    DLs.append(N_prop3[i][6])
    #    P_takeoffs.append(N_prop3[i][7])
    #
    #print('3 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #
    #for i in range(len(N_prop4)):
    #    RPMs.append(N_prop4[i][1])
    #    blade_amount.append(N_prop4[i][2])
    #    SPLs.append(N_prop4[i][3])
    #    Diameters.append(N_prop4[i][4])
    #    P_motors.append(N_prop4[i][5])
    #    DLs.append(N_prop4[i][6])
    #    P_takeoffs.append(N_prop4[i][7])
    #
    #print('4 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #
    #for i in range(len(N_prop5)):
    #    RPMs.append(N_prop5[i][1])
    #    blade_amount.append(N_prop5[i][2])
    #    SPLs.append(N_prop5[i][3])
    #    Diameters.append(N_prop5[i][4])
    #    P_motors.append(N_prop5[i][5])
    #    DLs.append(N_prop5[i][6])
    #    P_takeoffs.append(N_prop5[i][7])
    #
    #print('5 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    #
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #
    #for i in range(len(N_prop6)):
    #    RPMs.append(N_prop6[i][1])
    #    blade_amount.append(N_prop6[i][2])
    #    SPLs.append(N_prop6[i][3])
    #    Diameters.append(N_prop6[i][4])
    #    P_motors.append(N_prop6[i][5])
    #    DLs.append(N_prop6[i][6])
    #    P_takeoffs.append(N_prop6[i][7])
    #
    #print('6 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    #
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #
    #for i in range(len(N_prop7)):
    #    RPMs.append(N_prop7[i][1])
    #    blade_amount.append(N_prop7[i][2])
    #    SPLs.append(N_prop7[i][3])
    #    Diameters.append(N_prop7[i][4])
    #    P_motors.append(N_prop7[i][5])
    #    DLs.append(N_prop7[i][6])
    #    P_takeoffs.append(N_prop7[i][7])
    #
    #print('7 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    #
    #RPMs            = []
    #blade_amount    = [] 
    #Diameters       = [] 
    #SPLs            = []
    #P_motors        = []
    #DLs             = []
    #P_takeoffs      = []
    #
    #for i in range(len(N_prop8)):
    #    RPMs.append(N_prop8[i][1])
    #    blade_amount.append(N_prop8[i][2])
    #    SPLs.append(N_prop8[i][3])
    #    Diameters.append(N_prop8[i][4])
    #    P_motors.append(N_prop8[i][5])
    #    DLs.append(N_prop8[i][6])
    #    P_takeoffs.append(N_prop8[i][7])
    #
    #print('8 propellers Take-off/Hover \n'
    #      'Max RPM - Min RPM ','(',max(RPMs),' - ',min(RPMs),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount),' - ',min(blade_amount),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs),' - ',min(SPLs),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters),' - ',min(Diameters),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors),' - ',min(P_motors),') \n'
    #      'Max disk loading - Min disk loading ','(',max(DLs),' - ',min(DLs),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_takeoffs),' - ',min(P_takeoffs),') \n')
    
    #print(N_prop2)
    #print(parameter)             
    # =============================================================================
    # print('total possible combinations: ',len(parameter))
    # print('total possible combinations with 2 propellers: ',len(N_prop2))
    # print('total possible combinations with 3 propellers: ',len(N_prop3))
    # print('total possible combinations with 4 propellers: ',len(N_prop4))
    # print('total possible combinations with 5 propellers: ',len(N_prop5))
    # print('total possible combinations with 6 propellers: ',len(N_prop6))
    # print('total possible combinations with 7 propellers: ',len(N_prop7))
    # print('total possible combinations with 8 propellers: ',len(N_prop8))
    # =============================================================================
    
    # =============================================================================
    # |||||||||||||||||||||||||||||||CRUISE||||||||||||||||||||||||||||||||||||||||
    # =============================================================================
    
    PNL_A_cruise = noisedat.PNL_A_cruise -4   
    PNL_cruise = PNL_A_cruise + 14                                                                 #dB, requirement for aircraft at close range
    noy_cruise = 10**((PNL_cruise-26)/(33+1/3))                                                  #maximum value of noy to achieve the requirements for PNL A-weighted
    noy_range = np.arange(1,noy_cruise,0.07)
    acceptable_noys_c = [[50, 120],[63, 119],[80, 117],[100, 115],[125, 114],
                       [160, 113],[200, 111], [250, 110],[315,109],[400,108],
                       [500, 108],[630, 108],[800,108], [1000,108], [1250,106],
                       [1600,102],[2000,100],[2500,98],[3150,97],[4000,97],
                       [5000,98],[6300,99],[8000,102],[10000,105]]                    #[freq,SPL]
    
    N_blades_c    = np.arange(2,6,1)                                                   #constraint on number of blades. Min = 2, Max = 5
    parameter_c   = []                                                                 #To append in for loop next lines, [RPM,N_blades]                                                             #Speed of sound at sea level in m/s
    N_prop_c      = np.arange(1,2,1)                                                  #range of number of possible propellers
    
    
    
    k=1/(np.pi*A*ecc)
    WSstall=0.5*rhotrans*(Vstall**2)*CLmax
    Pcruise=(1/(Pset*mup))*(((CD0*0.5*rhoCruise*Vcruise**3)/(WSstall))+((WSstall)/(0.5*rhoCruise*Vcruise/k)))*(MTOW_c-0.25*Mpl)*9.81/1000
    #print(Pcruise)
    acc = 5.5  #m/s**2
    time_req_to_acc =  20
    P_req_to_acc = acc**2*2*MTOW*time_req_to_acc/1000
    
    N_prop1_c = []
    N_prop2_c = []
    N_prop3_c = []
    N_prop4_c = []
    N_prop5_c = []
    N_prop6N_c = []
    N_prop7_c = []
    N_prop8_c = []
    
    for i in range(len(N_blades_c)):
        for l in range(len(acceptable_noys_c)):
            RPM_c     = (acceptable_noys_c[l][0]*60)/N_blades_c[i]                          #RPMs
            #print(RPM_c)
            D_prop_c  = (M_tmax_c*60*sos_c)/(pi*RPM_c)                                       #b = diameter propeller
            #print(D_prop_c)
            if RPM_c > 1800 and RPM_c < 3300:   
                if D_prop_c < 5.5 and D_prop_c > 0.1:                                        #Define range of propeller diameter between 0.5 and 7.5m
                    for m in range(len(N_prop_c)):
                        P_mot_c               = e**(1/15.3*(acceptable_noys_c[l][1]-83.4+20*log(D_prop_c)-38.5*M_tmax_c+3*(N_blades_c[i]-2)-10*log(N_prop_c[m])+20*log(r_c)))    #Formula for SPL rewritten to get P_motor       
                        #print(P_mot_c)
                        P_mot_per_propeller_c = P_mot_c/N_prop_c[m]
                        #print(P_mot_c)
                        if Pcruise/N_prop_c[m] < 560*(N_motors500):                                    #Maximum technology motor now is 500kW so only take shit below this value
                            #DL_c          = MTOW_c*9.81/(pi*D_prop_c**2/4)/N_prop_c[m]    #Disk loading formula in N/m^2 per propeller
                            #WSstall=0.5*rhotrans*(Vstall**2)*CLmax
                            #Pcruise=(1/(Pset*mup))*(((CD0*0.5*rhoCruise*Vcruise**3)/(WSstall))+((WSstall)/(0.5*rhoCruise*Vcruise/k)))*(MTOW_c-0.25*Mpl)*9.81/1000
                            #P_req_hover_c = TW_to*sqrt(DL/(2*rhotrans*Ti**3))/FOM*MTOW*9.81/1000    #Power required for hover
                            #Vh_c          = sqrt(DL_c/(2*rhotrans))                          #Vh is induced velocity (or mass flow)
                            #P_to_c        = P_req_hover*((0.5*ROC_vtol/Vh)+sqrt(0.25*((ROC_vtol/Vh)**2+1)))
                            #print(Pcruise)
                            #if DL_c > 1000 and DL_c < 10000:  
                                         #Assumed disk loading between 100 and 1000 N/m^2 due to lower disk loading is more efficient for hover
                            #parameter_c.append(['#Propellers =',N_prop_c[m],'RPM =',RPM_c,'#Blades =',N_blades_c[i],'SPl =',acceptable_noys_c[l][1],'Diameter Prop =',D_prop_c,'P_mot =',P_mot_c,'Pcruise =',Pcruise])
                            
                            if Pcruise*1.25 < P_mot_c: 
                                if N_prop_c[m] > 0 and N_prop_c[m] < 2:
                                    Dm3 = round(D_prop_c,2)
                                    Dft3    = 3.2808399 * Dm3
                                    B_v    = N_blades_c[i]                  
                                    AF     = 150                    #Typical range is 100 - 150
                                    N_v    = RPM_c                  #prop speed (rpm)
                                    M = 0.22
                                    Pbr_kW_3 = Pcruise*1.25
                                    Pbr_hp_3 = 1.34102209 * Pbr_kW_3
                                    Kw     = 180
                                    W_prop_lbs3 = Kw*((Dft3/10)**2*(B_v/4)**0.7*(AF/1000)**0.75*(N_v*Dft3/20000)**0.5*(M+1)**0.5*(Pbr_hp_3/(10*Dft3**2))**0.12)
                                    W_prop_kg3  = 0.45359237 * W_prop_lbs3
                                    N_prop1_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise*1.25,2), round(W_prop_kg3*Pcruise*1.25,2), W_prop_kg3])
        #                        if N_prop_c[m] > 1 and N_prop_c[m] < 3:
        #                            N_prop2_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 2 and N_prop_c[m] < 4:
        #                            N_prop3_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 3 and N_prop_c[m] < 5:
        #                            N_prop4_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 4 and N_prop_c[m] < 6:
        #                            N_prop5_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 5 and N_prop_c[m] < 7:
        #                            N_prop6_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 6 and N_prop_c[m] < 8:
        #                            N_prop7_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #                        if N_prop_c[m] > 7 and N_prop_c[m] < 9:
        #                            N_prop8_c.append([N_prop_c[m],RPM_c,N_blades_c[i],acceptable_noys_c[l][1],round(D_prop_c,2),round(P_mot_c,2),round(Pcruise,2)])
        #
        #                        parameter_c.append(['#Propellers =',N_prop_c[m],'RPM =',RPM_c,'#Blades =',N_blades_c[i],'SPl =',acceptable_noys_c[l][1],'Diameter Prop =',D_prop_c,'P_mot =',P_mot_c,'Pcruise =',Pcruise])     #Appended in 'parameter' is [RPM,N_blades, SPL value, N_diameter, P_mot, Disk Loading, P_req, N_propellers] related to eachother
        #                        #print(P_mot)
    #print(len(parameter))
    
    #print('total possible combinations: ',len(parameter))
    #print('total possible combinations with 1 propellers: ',len(N_prop1_c))
    #print('total possible combinations with 2 propellers: ',len(N_prop2_c))
    #print('total possible combinations with 3 propellers: ',len(N_prop3_c))
    #print('total possible combinations with 4 propellers: ',len(N_prop4_c))
    #print('total possible combinations with 5 propellers: ',len(N_prop5_c))
    #print('total possible combinations with 6 propellers: ',len(N_prop6_c))
    #print('total possible combinations with 7 propellers: ',len(N_prop7_c))
    #print('total possible combinations with 8 propellers: ',len(N_prop8_c))
    
    
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop1_c)):
    #    RPMs_c.append(N_prop1_c[i][1])
    #    blade_amount_c.append(N_prop1_c[i][2])
    #    SPLs_c.append(N_prop1_c[i][3])
    #    Diameters_c.append(N_prop1_c[i][4])
    #    P_motors_c.append(N_prop1_c[i][5])
    #    P_cruise_c.append(N_prop1_c[i][6])
    #
    #print('1 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop2_c)):
    #    RPMs_c.append(N_prop2_c[i][1])
    #    blade_amount_c.append(N_prop2_c[i][2])
    #    SPLs_c.append(N_prop2_c[i][3])
    #    Diameters_c.append(N_prop2_c[i][4])
    #    P_motors_c.append(N_prop2_c[i][5])
    #    P_cruise_c.append(N_prop2_c[i][6])
    #
    #print('2 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop3_c)):
    #    RPMs_c.append(N_prop3_c[i][1])
    #    blade_amount_c.append(N_prop3_c[i][2])
    #    SPLs_c.append(N_prop3_c[i][3])
    #    Diameters_c.append(N_prop3_c[i][4])
    #    P_motors_c.append(N_prop3_c[i][5])
    #    P_cruise_c.append(N_prop3_c[i][6])
    #
    #print('3 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop4_c)):
    #    RPMs_c.append(N_prop4_c[i][1])
    #    blade_amount_c.append(N_prop4_c[i][2])
    #    SPLs_c.append(N_prop4_c[i][3])
    #    Diameters_c.append(N_prop4_c[i][4])
    #    P_motors_c.append(N_prop4_c[i][5])
    #    P_cruise_c.append(N_prop4_c[i][6])
    #
    #print('4 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop5_c)):
    #    RPMs_c.append(N_prop5_c[i][1])
    #    blade_amount_c.append(N_prop5_c[i][2])
    #    SPLs_c.append(N_prop5_c[i][3])
    #    Diameters_c.append(N_prop5_c[i][4])
    #    P_motors_c.append(N_prop5_c[i][5])
    #    P_cruise_c.append(N_prop5_c[i][6])
    #
    #print('5 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop6_c)):
    #    RPMs_c.append(N_prop6_c[i][1])
    #    blade_amount_c.append(N_prop6_c[i][2])
    #    SPLs_c.append(N_prop6_c[i][3])
    #    Diameters_c.append(N_prop6_c[i][4])
    #    P_motors_c.append(N_prop6_c[i][5])
    #    P_cruise_c.append(N_prop6_c[i][6])
    #
    #print('6 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop7_c)):
    #    RPMs_c.append(N_prop7_c[i][1])
    #    blade_amount_c.append(N_prop7_c[i][2])
    #    SPLs_c.append(N_prop7_c[i][3])
    #    Diameters_c.append(N_prop7_c[i][4])
    #    P_motors_c.append(N_prop7_c[i][5])
    #    P_cruise_c.append(N_prop7_c[i][6])
    #
    #print('7 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    #
    #
    #RPMs_c            = []
    #blade_amount_c    = [] 
    #Diameters_c       = [] 
    #SPLs_c            = []
    #P_motors_c        = []
    #P_cruise_c        = []
    #           
    #for i in range(len(N_prop8_c)):
    #    RPMs_c.append(N_prop8_c[i][1])
    #    blade_amount_c.append(N_prop8_c[i][2])
    #    SPLs_c.append(N_prop8_c[i][3])
    #    Diameters_c.append(N_prop8_c[i][4])
    #    P_motors_c.append(N_prop8_c[i][5])
    #    P_cruise_c.append(N_prop8_c[i][6])
    #
    #print('8 propellers Cruise \n'
    #      'Max RPM - Min RPM ','(',max(RPMs_c),' - ',min(RPMs_c),') \n' 
    #      'Max blades - Min blades ','(',max(blade_amount_c),' - ',min(blade_amount_c),') \n'
    #      'Max SPL - Min SPL ','(',max(SPLs_c),' - ',min(SPLs_c),') \n'
    #      'Max diameter - Min diameters ','(',max(Diameters_c),' - ',min(Diameters_c),') \n'
    #      'Max P_motors - Min P_motors ','(',max(P_motors_c),' - ',min(P_motors_c),') \n'
    #      'Max P_takeoff - Min P_takeoff ','(',max(P_cruise_c),' - ',min(P_cruise_c),') \n')
    
    listdiameter = []
    
    for i in range(len(N_prop4)):
        listdiameter.append(N_prop4[i][4])
    
    combinations = []
    Diameter = min(listdiameter)
    for i in range(len(N_prop4)):
        if Diameter == N_prop4[i][4]:
            combinations.append(N_prop4[i])
                  
    #print (combinations)
    
    listdiameter_c = []
    
    for i in range(len(N_prop1_c)):
        listdiameter_c.append(N_prop1_c[i][7])
      
    combinations_c = []
    PowerWeight = min(listdiameter_c)
    for i in range(len(N_prop1_c)):
        if PowerWeight == N_prop1_c[i][7]:
            combinations_c.append(N_prop1_c[i])
    #print(combinations_c)
    
            
    
    Area = 4*(Diameter/2)**2*pi
    Area_1 = Area*0.6
    Area_2 = Area*0.4
    
    Dm1     = sqrt(Area_1/(pi*2))*2                 #single propeller system diameter
    Dft1    = 3.2808399 * Dm1         #3.2808399 multiplication is to go from m to ft.
    Dm2     = sqrt(Area_2/(pi*2))*2                   #single propeller system diameter
    Dft2    = 3.2808399 * Dm2         #3.2808399 multiplication is to go from m to ft.
    #Dft3    = 3.2808399 * Dm3 
    
    B      = combinations[0][2]                      #number of blades
    #B_v    = 5                  
    AF     = 150                    #Typical range is 100 - 150
    N      = combinations[0][1]                   #prop speed (rpm)
    #N_v    = 3750
    M      = 0.22                   #airplane mach speed
    Pbr_kW_1 = combinations[0][-1]  
    Pbr_kW_2 = combinations[0][-1]                  #shaft power in kW
    #Pbr_kW_3 = 236                
    Pbr_hp_1 = 1.34102209 * Pbr_kW_1
    Pbr_hp_2 = 1.34102209 * Pbr_kW_2    #1.34102209 multiplication is to go from kW to horse power.
    #Pbr_hp_3 = 1.34102209 * Pbr_kW_3
    Kw     = 180                    #range of constant is 160-180 for composite propellers
    
    W_prop_lbs1 = Kw*((Dft1/10)**2*(B/4)**0.7*(AF/1000)**0.75*(N*Dft1/20000)**0.5*(M+1)**0.5*(Pbr_hp_1/(10*Dft1**2))**0.12)
    W_prop_lbs2 = Kw*((Dft2/10)**2*(B/4)**0.7*(AF/1000)**0.75*(N*Dft2/20000)**0.5*(M+1)**0.5*(Pbr_hp_2/(10*Dft2**2))**0.12)
    #W_prop_lbs3 = Kw*((Dft3/10)**2*(B_v/4)**0.7*(AF/1000)**0.75*(N_v*Dft3/20000)**0.5*(M+1)**0.5*(Pbr_hp_3/(10*Dft3**2))**0.12)
    
    W_prop_kg1  = 0.45359237 * W_prop_lbs1
    W_prop_kg2  = 0.45359237 * W_prop_lbs2
    #W_prop_kg3  = 0.45359237 * W_prop_lbs3
    
    D_prop_big = Dm1 #This is the diameter of the big vertical propeller
    DL_total = (combinations[0][6]) #This is the total disk Loading
    D_prop_small = (Dm2)#This is the diameter of the small vertical propeller
    D_prop_ff = (combinations_c[0][4]) #This is the diameter of the ff horizontal propeller
    
    #print(W_prop_kg1) #THIS IS THE WEIGHT OF THE BIG VERTICAL PROPELLER (ONE)
    #print(W_prop_kg2) #THIS IS THE WEIGHT OF THE SMALL VERTICAL PROPELLER (ONE)
    #print(combinations_c[0][8]) #THIS IS THE WEIGHT OF THE FF HORIZONTAL PROPELLER (ONE)
    
    WProp_big = W_prop_kg1 + 12*2 + 72
    WProp_small = W_prop_kg2 + 12*2 + 72
    Wprop_ff = combinations_c[0][8] + 12*4*N_motors500 + 135*N_motors500
    RPM_vtol = (combinations[0][1])
    RPM_cruise = (combinations_c[0][1])
    N_blades_vtol = (combinations[0][2])
    N_blades_cruise = (combinations_c[0][2])
#    print(N_blades_vtol, N_blades_cruise)
#    print(RPM_vtol,RPM_cruise)
    if test==None:    
        return(WProp_big, WProp_small,Wprop_ff,DL_total,D_prop_ff,D_prop_big,D_prop_small, Pbr_kW_1*1000, combinations_c[0][-3]*1000)
    else:
        return (D_prop_small, D_prop_big, N_blades_vtol, N_blades_cruise, RPM_vtol, RPM_cruise,  Pbr_kW_1)




    
    