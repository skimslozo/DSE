import numpy as np
from readPars import load_input, plot_graph
#aero = load_input(sheet = 'Aero')


def reynold(rho, V_freestream, L, mu):
    number = rho*V_freestream*L/mu
    return number

def design(V, rho, W, S1, S2, S3, sweep1, sweep2, sweep3, factor_tail, Cl_1, Cl_2):
    '''Wing'''
    q = 0.5 * rho * V**2
    CL_design = (1 + factor_tail)/q*(W/S) #cl design calculation; assumption; decrease of hydrogen almost no influence --> not average 
    
    '''Airfoil'''
    #airfoil centerbody theoretical
    Veff1 = V*np.cos(sweep1)
    qeff1 = .5*rho*Veff1**2
    Cl1 = CL_design*q/qeff1
    
    #airfoil transition theoretical
    Veff2 = V*np.cos(sweep1)
    qeff2 = .5*rho*Veff2**2
    Cl2 = CL_design*q/qeff2
    
    #airfoil outer theoretical
    Veff3 = V*np.cos(sweep1)
    qeff3 = .5*rho*Veff3**2
    Cl3 = CL_design*q/qeff3
    
    #scaling to achieve enough L over the whole wing
    SCl = S1*Cl1 + S2*Cl2 + S3*Cl3
    Cl_3 = (SCl - S1*Cl_1 - S2*Cl_2)/S3
    
    return Cl_1, Cl_2, Cl_3

mission = load_input(sheet = 'Mission')
sweep_wing =  0#rad ---------------------------------------------
tailfactor = .1
V_freestream = mission.Vcruise
c1_avg = 9.0
c2_avg = 8. 
c3_avg = 12.8/2
#make distinction between different airfoils

number = 1.225*250/3.6*9/(1.81*10**-5)/1000000
#V = mission.Vcruise, rho=mission.rhoCruise, 
#design(mission.Vcruise, mission.rhoCruise, ..., ..., tailfactor)