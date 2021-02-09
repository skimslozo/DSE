# -*- coding: utf-8 -*-
"""
Created on Mon May 06 16:27:34 2019

@author: Miks
"""

import matplotlib.pyplot as plt
import numpy as np
from DesIter_H2_Ray import IterNew


Npropsrange=np.array([4.])#np.arange(2, 6, 2)
Nbladesrange=np.array([6.])
Arange=np.arange(1, 3.5, 0.5)
DLrange=np.array([976*4.])
Vcruiserange=np.arange(250, 320, 20)                
wmin=1e20
niter=60
iters=np.linspace(1, niter, niter)
MTOWi=3170.
weights=[]

for Nprops in Npropsrange:
    for Nblades in Nbladesrange:
        for A in Arange:
            for DL in DLrange:
                for Vcruise in Vcruiserange:
                    pars=[Nprops, Nblades, A, DL, Vcruise]
                    print(pars)
                    for i in iters:
                        if i==1:
                            ans=IterNew(MTOWi=MTOWi, pars=pars, mode='optimize')
                        else:
                            ans=IterNew(MTOWi=ans[0], pars=pars, mode='optimize')
                        if ans[0]>ans[1]:
                            break
                    if ans[0]<wmin and ans[0]<=ans[1] and ans[0]>0 and ans[5]>0.04:
                        wmin=ans[0]
                        bestpars=pars
                        Surface=ans[2]
                        Aspect=ans[3]
                        Wingspan=ans[4]
                        propspace=ans[5]
                        Dprop=ans[6]
                        ROCvtol=ans[8]
                        DLvtol=ans[9]
                        Vcmin=ans[10]
                        Ratio=ans[11]
print('Best design:')
print('Number of propellers: {}'.format(bestpars[0]))
print('Number of blades: {}'.format(bestpars[1]))
print('Aspect ratio: {}'.format(bestpars[2]))
print('Wing surface area: {} m^2'.format(Surface))
print('Maximum take-off weight: {} kg'.format(wmin))
print('Propeller spacing: {} m'.format(propspace))
print('Wingspan: {} m'.format(Wingspan))
print('Propeller diameter: {} m'.format(Dprop))
print('ROCvtol: {} m/s'.format(ROCvtol))
print('Disk loading: {} N/m2'.format(DLvtol))
print('Cruise speed: {} km/h'.format(Vcmin*3.6))
print('Pcruise/Pto: {} '.format(Ratio))