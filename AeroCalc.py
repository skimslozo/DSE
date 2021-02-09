# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 23:35:30 2019

@author: timvd
"""
import numpy as np

def getCD0(Swing, Stail): #lift and drag estimation page 42
    cdc_wing = 0.007
    cdc_fuselage = .11
    cdc_tail = .008
    cdmiscfactor = 1.15
    L1 = 2
    L2 = 6
    L3 = 3
    D = 1.5
    Sfuselage = np.pi*D/4*(1/(3*L1**2)*(((4*L1**2 + .25*D**2)**1.5 - 1/8*D**3)) - D + 4*L2 + 2*(L3**2 + .25*D**2)**(0.5))/2+ 1.5* (L1 + L2 + L3) #only half round fuselage, also partly wing 
    Sref = (Sfuselage + Swing + Stail)*2.9
    CD0 = 1/Sref*(cdc_wing*Swing*2*1.07 + (cdc_fuselage + cdc_wing)/2*Sfuselage + cdc_tail*Stail)*cdmiscfactor
    return CD0*1.1

#cmac= -.04 #given based on airfoil analysis, doesn't change
#twist = 0 #this has been based on an article and difficulty to implement this based on the given design configuration, ask if explanation is necessary
#AR = 
#cla = (.8+.5)/np.radians(6-(-6)) #if desired ask me for explanation
#M = 
#sweephalf = #see below in comments for def for it
#sweepc4 = 
#cl0 =0.181 #airfoil data
#cma = 0 #based on data
#
#"""slope of CL curve calculation"""
#def CLA(AR, cla, M, sweephalf): #(aspect ratio, cl_alpha airfoil, machnumber, sweep of half chord) #all in radians, datcom
#    beta = np.sqrt(1-M**2)
#    k = cla/(2*np.pi)
#    CL_a = (2*np.pi * AR)/(2 + np.sqrt((AR*beta/k)**2 * (1+ np.tan(sweephalf)**2 / beta**2) + 4)) #346 method
#    return CL_a
#CL_a = CLA(AR, cla, M, sweephalf) #verified with example on 349 snorri #radians
#
#"""Clmax calculation; this one is harder to implement as manual work is required; ask me how it works and what to input"""
#def CLMAX(twist, clmax, sweepc4, CL_a, AR): #(overall twist of wing, clmax of airfoil, sweep c/4, CLa of wing; function above, aspect ratio)  philips and alley method 
#    ktaper1 = float(input("What is ktaper1?")) #page 13
#    ktaper2 = float(input("What is ktaper2?"))
#    CLtwist0_clmaxtwist0 = float(input("What is CLtwist0_clmaxtwist0?"))#philips and alleys page 7
#    kltaper = 1 + ktaper1*sweepc4 - ktaper2*sweepc4**1.2
#    #kltwist = (1- (CL_clmax)/(CLtwist0_clmaxtwist0))/CL_a*twist/clmax
#    kltwist = float(input("What is kltwist? CLalpha*twist/clmax")) #page 10
#    kls = 1 + (0.0042*AR - 0.068)*(1+2.3*CL_a*twist/clmax)
#    CLmax = CLtwist0_clmaxtwist0 *kls*kltaper*(clmax - kltwist*CL_a*twist)
#    return CLmax
#CLmax = CLMAX(twist, clmax, SWEEPC4(sweepLE, rootchord, span, taper), CL_a, AR)
#
#"""3d characteristics of 2d airfoil"""
#def characteristics3d(cl0, CLalfa, clalfa, cmalfa): #gives 3d characteristics of airfoil data #(cl0 airfoil, CLa wing, cla 2d airfoil, cma airfoil)                  #all verified example 349 snorri
#    alfa0 = -cl0/clalfa #radians
#    CL0   = -alfa0*CLalfa #close to datcom method snorri
#    Cmalfa = CLalfa*(cmalfa/clalfa)
#    return [alfa0, CL0, Cmalfa]
#alfa0, CL0, Cmalfa = characteristics3d(cl0, CL_a, cla, cmalfa)[0], characteristics3d(cl0, CL_a, cla, cmalfa)[1], characteristics3d(cl0, CL_a, cla, cmalfa)[2]
#
#
#""" Very general CD0 calculation; based on adsee slides"""
##def CD0(cdc_wing, Swing, cdc_fuselage, cdc_tail, Stail, Sref, cdmiscfactor, D, L1, L2, L3 ): #lift and drag estimation page 42
##    cdcwing = 0.007
##    cdcfuselage = .11
##    cdctail = .008
##    L1 = 2 
##    L2 = 6
##    L3 = 3
##    cdmiscfactor = 1.15
##    Sfuselage = np.pi*D/4*(1/(3*L1**2)*(((4*L1**2 + .25*D**2)**1.5 - 1/8*D**3)) - D + 4*L2 + 2*(L3**2 + .25*D**2)**(0.5))/2+ 1.5* (L1 + L2 + L3) #only half round fuselage, also partly wing 
##    CD0 = 1/Sref*(cdc_wing*Swing*2*1.07 + (cdc_fuselage + cdc_wing)/2*Sfuselage + cdc_tail*tailarea)*cdmiscfactor
##    return CD0 
#
#
#"""CD calc"""
#def CD(CD0, CL, AR, sweepLE): #snorri
#    #e = 4.61*(1 - 0.045*AR**.68)*(np.cos(sweepLE)**0.15) - 3.1 #MohammadH.Sadraey-AircraftDesignASystemsEngineeringApproach-JohnWileyamp_Sons2012 pp. 212--> this one doesn't work for this taper and sweep
#    e = 1.78*(1 - 0.045*AR**.68) - 0.64
#    CD = CD0 + (CL**2) / (np.pi * AR * e) #p.693 snorri
#    return [CD, e]
## =============================================================================
## def sweephalfway(span, rootchord, taper, sweepLE): 
##     tipchord = taper * rootchord
##     vertical_LE = span/2*  np.tan(sweepLE)
##     vertical_TE = rootchord - vertical_LE - tipchord
##     sweepTE = np.arctan(vertical_TE/(span/2))
##     sweephalf = (sweepLE + sweepTE)/2
##     return [sweephalf, sweepTE]
## =============================================================================

