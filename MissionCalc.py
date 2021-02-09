# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:50:41 2019

@author: Miks
"""

import numpy as np
import matplotlib.pyplot as plt
from readPars import load_input

#MTOWi=3353.9570706384593
#pars=[4., 2., 5.0, 1300., 250.]

missiondat=load_input(sheet='Mission')
acdat=load_input(sheet='Aircraft parameters')
energdat=load_input(sheet='Energy storage')
miscdat=load_input(sheet='Misc')

def getClassII(MTOWi=acdat.MTOW, DL=acdat.DL, CD0in=0.01, pars=None, missiondat=missiondat, acdat=acdat, energdat=energdat, miscdat=miscdat, retrn=None, mode='read'):
    '''
       retrn - Variabl name as a string that you want to output in read mode
       mode - 'read' or 'optimize' - optimize only used if new parameters have to be re-iterated with new design changes
    '''
    showdiag=False
    #========================================================================================
    #Known parameter values (Class I)
    #========================================================================================
    Espec=energdat.Espec
    WmotPmax=miscdat.WmotPmax
    if mode=='read':
       Nprops=acdat.Nprops
       Nblades=acdat.Nblades
       A=acdat.A
       Vcruise=missiondat.Vcruise
       MTOWi=acdat.MTOW
    elif mode=='optimize':
        retrn=None
        if pars==None:
            raise Exception('Optimization parameters not specified')
        Nprops=pars[0]
        Nblades=pars[1]
        A=pars[2]    
        DL=pars[3]
        Vcruise=pars[4]/3.6
    else:
        raise Exception('Please specify the usage of the function in mode input parameter')
    
    rhoSL=missiondat.rhoSL #Density at sea level
    rhoTrans=missiondat.rhoTrans #Density at transition altitude of hT=500m
    rhoCruise=missiondat.rhoCruise #Density at cruise altitude of hC=3000ft
    bmax=acdat.bmax
    #========================================================================================
    #Assumed parameter values (Class I)
    #========================================================================================
    
    CLmax=acdat.CLmax
#    CD0=acdat.CD0
    CD0=CD0in
    Vstall=0.5*Vcruise
    Pset=acdat.Pset
    mup=acdat.mup
    e=acdat.e
    k=1/(np.pi*A*e)
    
    ROCff=missiondat.ROCff
    CLROC=np.sqrt(3*CD0*1/k)
    CDROC=4*CD0
    
    M=acdat.M 
    TWto=acdat.TWto
    TWland=acdat.TWland
    Ti=acdat.Ti
    #========================================================================================
    #Requirement value generation
    #========================================================================================
    WSgen=np.arange(1, 3000)
    WSstall=(0.5*rhoTrans*(Vstall**2)*CLmax/1.3)
    PWcruise=(1/(Pset*mup))*(((CD0*0.5*rhoCruise*Vcruise**3)/(WSgen))+((WSgen)/(0.5*rhoCruise*Vcruise/k)))
    PWroc=(1/mup)*(ROCff+((np.sqrt(WSgen*(2/rhoCruise)))/((CLROC**1.5)/CDROC)))
    PWrod=0.5*rhoCruise*(Vcruise**3)*(1./WSgen)*(CD0+k*WSgen/(0.5*rhoCruise*Vcruise**2))
    PWhover=TWto*np.sqrt((DL)/(2*rhoTrans*(Ti**3)))*(1./M)
    PWland=TWland*np.sqrt((DL)/(2*rhoTrans*(Ti**3)))*(1./M) # Not sure if the rate of descent foruma applicable here, so we're just going for the rate of descent formula, with T/W=1
    
    Smax=(bmax**2)/A

    Kmaterial=miscdat.Kmaterial
    Hcruise=missiondat.Hcruise
    R=missiondat.R
    Mpl=missiondat.Mpl
    HHVH2=energdat.HHVH2 #Higher heating value of H2, J/kg
    
    N_mot_hor=1
    N_prop_hor=N_mot_hor
    #========================================================================================
    #Assumed parameter values (Class II)
    #========================================================================================
    Kprop=miscdat.Kprop
    Htrans=missiondat.Htrans
    DoD=energdat.DoD
    MFstruct=acdat.MFstruct
    Msubsys=acdat.Msubsys
    Mavion=acdat.Mavion
    mubat=energdat.mubat
    mu_comp=energdat.mu_comp #Energy loss due to pump/valve/compressor losses
    mupem=energdat.mupem # PEM stack efficiency
    mutank=energdat.mutank #Tank efficiency estimate
    tankfactor=energdat.tankfactor #Some factor to take mass increase due to some tank component into account
    mfact_comp=energdat.mfact_comp #30% to take into account the mass of the fuel cell + battery components   
    
    ROCvtol=missiondat.ROCvtol
    ROCfactor=missiondat.ROCfactor
    Wcabin=1.56
    #========================================================================================
    #Estimation of the parameters
    #========================================================================================
    Phover=PWhover*MTOWi*9.81
    S=MTOWi*9.81/WSstall
    b=np.sqrt(A*S)
    if b>bmax:
        raise Exception('Wingspan exceeds maximum allowable wingspan')
    ROCh=ROCfactor*ROCff
    Wmax=WSstall*Smax/9.81
    vh=np.sqrt(DL/(2*rhoTrans))
    xyz=((0.5*ROCvtol/vh)+np.sqrt(0.25*((ROCvtol/vh)**2)+1))
    Pto=((0.5*ROCvtol/vh)+np.sqrt(0.25*((ROCvtol/vh)**2)+1))*Phover
    Pclimb=(1/mup)*(ROCh+((np.sqrt(WSstall*(2./rhoCruise)))/((CLROC**1.5)/CDROC)))*MTOWi*9.81
    Pcruise=(1/(Pset*mup))*(((CD0*0.5*rhoCruise*Vcruise**3)/(WSstall))+((WSstall)/(0.5*rhoCruise*Vcruise/k)))*(MTOWi-0.25*Mpl)*9.81
    Pdesc=(0.5*rhoCruise*(Vcruise**3)*(1./WSstall)*(CD0+k*WSstall/(0.5*rhoCruise*Vcruise**2)))*MTOWi*9.81 #Assumption - descent in forward flight is performed at cruise velocity
    if -ROCvtol/vh<=-2:
        Pland=((0.5*-ROCvtol/vh)-np.sqrt(0.25*((-ROCvtol/vh)**2)-1))*Phover
    else:
        Pland=((-ROCvtol/vh)+0.974-1.125*(-ROCvtol/vh)-1.372*((-ROCvtol/vh)**2)*-1.718*((-ROCvtol/vh)**3)-0.655*((-ROCvtol/vh)**4))*Phover
    VROC=(Pclimb/(0.5*rhoCruise*S*CDROC))**(1/3.)
    
    Tvtol=Htrans/ROCvtol
    Tclimb=(Hcruise-Htrans)/ROCh
    Tcruise=((R/3.)-VROC*np.cos(np.arcsin(ROCh/VROC))*Tclimb-Vcruise*np.cos(np.arcsin(ROCh/VROC))*Tclimb)/Vcruise #Insert VROC
    Tsegment=Tvtol*2+Tclimb*2+Tcruise
    Tmission=Tsegment*3
    
    CLcruise=(MTOWi*9.81)/(0.5*rhoCruise*S*Vcruise**2)
    CDcruise=CD0+k*CLcruise**2
    LD=CLcruise/CDcruise
    
    return locals()[retrn]
