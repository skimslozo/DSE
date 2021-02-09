# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:57:40 2019

@author: Miks
"""

import numpy as np
from readPars import load_input, plot_graph
import matplotlib.pyplot as plt
from DesIter_H2_Ray import IterNew
import seaborn as sns

class Geometry:
    
    def __init__(self):
        self.acdat=load_input(sheet='Aircraft parameters')
        self.taildat=load_input(sheet='Tail parameters')
        self.geodat=load_input(sheet='Dimensions and Positions')
        self.missiondat=load_input(sheet='Mission')
        self.b=IterNew(mode='read', retrn='b')
        self.bh=np.sqrt(self.taildat.Sh*self.taildat.Ah)
        self.S=IterNew(mode='read', retrn='S')
        self.Sh=self.taildat.Sh
        self.A=IterNew(mode='read', retrn='A')
        self.Ah=self.taildat.Ah
        self.xCr=self.geodat.xCr
        self.taper=0.2*(2-np.radians(self.acdat.lambdac4))
        self.taperh=self.taildat.taperh
        self.rootchord=(2*self.S)/((1+self.taper)*self.b)
        self.tipchord=self.taper*self.rootchord
        self.rch=(2*self.Sh)/((1+self.taperh)*self.bh)
        self.tch=self.taperh*self.rch
#FIX THIS!!!!!!!!!!111!        self.lambdac4h=np.arctan(self.bh/(self.rch-self.tch))
        self.lambdac4h=15.
        self.lambdaLE=np.arctan(np.tan(np.radians(self.acdat.lambdac4))+((1-self.taper)/(self.A*(1+self.taper))))
        self.lambdaLEh=np.arctan(np.tan(np.radians(self.lambdac4h))+((1-self.taperh)/(self.taildat.Ah*(1+self.taperh))))
        self.MGC=(2./3.)*self.rootchord*((1+self.taper+self.taper**2)/(1+self.taper))
        self.MGCh=(2./3.)*self.rch*((1+self.taperh+self.taperh**2)/(1+self.taperh))
        self.yMGC=(self.b/6)*((1+2*self.taper)/(1+self.taper))
        self.yMGCh=(self.bh/6)*((1+2*self.taperh)/(1+self.taperh))
        self.xMGC=self.yMGC*np.tan(self.lambdaLE)+self.xCr
        self.Nprops=int(IterNew(mode='read', retrn='Nprops'))
        self.Npropshor=int(IterNew(mode='read', retrn='N_mot_hor'))
        self.Dprop=IterNew(mode='read', retrn='Dprop')
        self.Dprophor=IterNew(mode='read', retrn='Dprop_hor')
        self.enersys=IterNew(mode='read', retrn='Menersys')
        self.fusclear=0.0
        self.LEclear=0
        self.interpropclear=0.0
        self.plotCGs=True
        self.Pfractions=[0.15, 0.35, 0.15, 0.35]
        self.Mprops=[self.Dprop/self.Nprops for i in range(self.Nprops)]
        self.show=True
        if sum(self.Pfractions)!=1:
            raise Exception('VTOL power distribution does not sum up to 1!')
        self.propinfo=[]
        self.Mbigprop=51.31
        self.Msmallprop=32.44
        self.Dbigprop=3.824
        self.Dsmallprop=3.122
        self.Dhorprop=1.38
        
    def getPropDiameter(self, Pfraction):
        self.Sprop=(Pfraction*IterNew(mode='read', retrn='MTOWi')*9.81)/(IterNew(mode='read', retrn='DL'))
        self.Dprop=np.sqrt(4*self.Sprop/np.pi)
        return self.Dprop
    
    def getPositions(self):
        self.yproplead=self.geodat.widthfus/2+self.getPropDiameter(self.Pfractions[0])/2+self.fusclear
        self.yproptrail=self.geodat.widthfus/2+self.getPropDiameter(self.Pfractions[1])/2+self.fusclear
        self.xleadprop=self.xCr+self.yproplead*np.tan(self.lambdaLE)+self.getPropDiameter(self.Pfractions[0])/2+self.LEclear
        self.proppos=[[self.xleadprop, self.yproplead], 
                      [self.xleadprop+self.getPropDiameter(self.Pfractions[1])/2+self.getPropDiameter(self.Pfractions[0])/2+self.interpropclear, self.yproptrail],
                      [self.xleadprop, -(self.yproplead)], 
                      [self.xleadprop+self.getPropDiameter(self.Pfractions[1])/2+self.getPropDiameter(self.Pfractions[0])/2+self.interpropclear, -(self.yproptrail)]]
        for i in range(self.Nprops):
            self.propinfo.append([self.proppos[i][0], self.proppos[i][1], self.Pfractions[i]])
        
    def addDiameters(self):
        for i in range(self.Nprops):
            self.propinfo[i].append(self.getPropDiameter(self.propinfo[i][2]))
            
    def getWingGeometry(self):
        self.xwing=np.array([0, self.b/2*np.tan(self.lambdaLE)])
        self.ywing=np.array([0, self.b/2])
        self.xcabin=np.linspace(0, self.geodat.lfus)
        self.xwing=np.append(self.xwing, [self.xwing[-1]+self.tipchord, self.rootchord])
        self.ywing=np.append(self.ywing, [self.b/2, 0])
        
        self.xwing=np.append(self.xwing, self.xwing)
        self.xwing=self.xwing+self.xCr
        self.ywing=np.append(self.ywing, -self.ywing)
        
    def getTailGeometry(self):
        self.xtail=np.array([0, self.bh/2*np.tan(self.lambdaLEh)])
        self.ytail=np.array([0, self.bh/2])
        self.xtail=np.append(self.xtail, [self.xtail[-1]+self.tch, self.rch])
        self.ytail=np.append(self.ytail, [self.bh/2, 0])
        
        self.xtail=np.append(self.xtail, self.xtail)
        self.xtail=self.xtail+self.xCrh
        self.ytail=np.append(self.ytail, -self.ytail)

                
    def getPropMass(self, Dprop, Pfraction):
        Mprop=IterNew(mode='read', retrn='Kmaterial')*IterNew(mode='read', retrn='Kprop')*(IterNew(mode='read', retrn='Nblades')**(0.391))*(Dprop*Pfraction*IterNew(mode='read', retrn='Pto')/(1000))**0.782
        MESC=0.7383e-04*(IterNew(mode='read', retrn='Pto')*Pfraction)**(0.8854)
        Mpropulsion=1.1*((IterNew(mode='read', retrn='WmotPmax')*((IterNew(mode='read', retrn='Pto')*Pfraction)/(9.81))+MESC)+Mprop)
        return Mpropulsion
    
    def getFFPropMass(self, Dprop):
        Mprop=self.Npropshor*IterNew(mode='read', retrn='Kmaterial')*IterNew(mode='read', retrn='Kprop')*(IterNew(mode='read', retrn='Nblades')**(0.391))*(Dprop*IterNew(mode='read', retrn='Pcruise')/(1000*self.Npropshor))**0.782
        MESC=0.7383e-04*(IterNew(mode='read', retrn='Pcruise')/self.Npropshor)**(0.8854)
        Mpropulsion=1.1*(self.Npropshor*(IterNew(mode='read', retrn='WmotPmax')*((IterNew(mode='read', retrn='Pcruise')/self.Npropshor)/(9.81))+MESC)+Mprop)
        return Mpropulsion
        
    def getCabinGeometry(self):
        self.xcab=np.array([(self.geodat.widthfus/2)*(1-np.cos(i)) for i in np.linspace(0, np.pi/2, 70)])
        self.ycab=np.array([(self.geodat.widthfus/2)*(np.sin(i)) for i in np.linspace(0, np.pi/2, 70)])

        self.xcab=np.append(self.xcab, [self.geodat.lfus, self.geodat.lfus])
        self.ycab=np.append(self.ycab, [self.geodat.widthfus/2, 0])
        
        self.xcab=np.append(self.xcab, self.xcab)
        self.ycab=np.append(self.ycab, -self.ycab)
        
    def getVTOLPropCG(self):
        self.Mvtolprop=0
        self.propx=0
        self.propdraw=[]
        for i in range(self.Nprops):
            masspr=self.getPropMass(self.propinfo[i][3], self.propinfo[i][2])
            self.propx+=self.propinfo[i][0]*masspr
            self.Mvtolprop+=masspr
            self.propdraw.append(self.getCircle(self.propinfo[i][3], [self.propinfo[i][0], self.propinfo[i][1]]))
        
    def getChordLength(self, y):
        return self.rootchord-((self.rootchord-self.tipchord)/(self.b/2))*abs(y)
    
    def getTailChord(self, y):
        return self.rch-((self.rch-self.tch)/(self.bh/2))*abs(y)

        
    def getComponentWeights(self):
        self.Swingw=self.S-2*(0.5*(self.rootchord+self.getChordLength(self.geodat.widthfus/2))*self.geodat.widthfus/2)
        self.Mwing=0.453592*(0.0051*((1.5*self.acdat.limitload*IterNew(mode='read', retrn='MTOWnew')/0.453592)**0.557)*((self.Swingw/(0.3048**2))**0.649)*(self.A**0.5)*((self.acdat.tcratio)**(-0.4))*((1+self.taper)**0.1)*(1/np.cos(self.acdat.lambdac4))*(self.geodat.Scsw**0.1))
        self.Kws=0.75*((1+2*self.taper)/(1+self.taper))*((self.b*np.tan(self.acdat.lambdac4)/self.geodat.lfus))
        self.Mfus=0.453592*(0.328*1.06*1.12*((1.5*self.acdat.limitload*IterNew(mode='read', retrn='MTOWnew')/0.453592)**0.5)*((self.geodat.lfus/0.3048)**0.25)*(((self.geodat.lfus*np.pi*(self.geodat.widthfus/2)**2)/(0.3048**2))**0.302)*((1+self.Kws)**0.04)*(IterNew(mode='read', retrn='LD')**0.1))
        self.Mhortail=0.453592*(0.0379*((1+self.geodat.widthfus/self.b)**(-0.25))*((IterNew(mode='read', retrn='MTOWnew')/0.453592)**0.639)*((1.5*self.acdat.limitload)**0.1)*((self.taildat.Sh/(0.3048**2))**0.75)*(1./(self.geodat.lh/0.3048))*((0.3*self.geodat.lh/0.3048)**0.704)*(1/np.cos(np.radians(self.lambdac4h)))*(self.taildat.Ah**0.166)*((1+self.taildat.SeSh)**0.1))
        self.Mverttail=1.7*self.Mhortail
        self.Menersys=IterNew(mode='read', retrn='Menersys')
        
        
    def getCGs(self):
        self.xcgwing=self.xCr+0.35*(self.b/2)*np.tan(self.lambdaLE)+(self.rootchord-0.35*(self.rootchord-self.tipchord))*0.585
        self.xcgfus=self.geodat.xcgfus*self.geodat.lfus
        self.propxcg=self.propx/self.Mvtolprop
        self.empxcg=self.xcgh
        self.xcgpl=(self.missiondat.Mpax*self.geodat.lfus*(2*self.geodat.Xcgpilot+self.geodat.Xcgseat3+self.geodat.Xcgseat4)+self.missiondat.Mpatient*self.geodat.Xcgpatient*self.geodat.lfus+self.geodat.Xcgeq*self.geodat.lfus*self.missiondat.Meq)/(4*self.missiondat.Mpax+self.missiondat.Mpatient+self.missiondat.Meq)
        self.mpl=4*self.missiondat.Mpax+self.missiondat.Mpatient+self.missiondat.Meq
#        self.xcgpl=(self.geodat.Xcgeq*self.geodat.lfus*self.missiondat.Meq)/(self.missiondat.Meq)
#        self.mpl=self.missiondat.Meq

        self.propFFmass=self.getFFPropMass(IterNew(mode='read', retrn='Dprop_hor'))
        self.xcgenersys=self.geodat.lfus*self.geodat.Xcgenersys
        
        self.totalxcg=(self.xcgwing*self.Mwing+self.xcgfus*self.Mfus+self.propxcg*self.Mvtolprop+self.empxcg*(self.Mhortail+self.Mverttail)+self.xffprop*self.propFFmass+self.xcgpl*self.mpl+self.xcgenersys*self.Menersys)/(self.Menersys+self.mpl+self.propFFmass+self.Mvtolprop+self.Mfus+self.Mwing+self.Mhortail+self.Mverttail)
        self.totalmass=self.Menersys+self.mpl+self.propFFmass+self.Mvtolprop+self.Mfus+self.Mwing+self.Mhortail+self.Mverttail
        
        self.xcgoew=(self.xcgwing*self.Mwing+self.xcgfus*self.Mfus+self.propxcg*self.Mvtolprop+self.empxcg*(self.Mhortail+self.Mverttail)+self.xffprop*self.propFFmass+self.geodat.Xcgeq*self.geodat.lfus*self.missiondat.Meq+self.xcgenersys*self.Menersys+self.geodat.Xcgpilot*self.geodat.lfus*self.missiondat.Mpax)/(self.missiondat.Mpax+self.missiondat.Meq+self.Menersys+self.propFFmass+self.Mvtolprop+self.Mfus+self.Mwing+self.Mhortail+self.Mverttail)
        self.Moew=self.missiondat.Mpax+self.missiondat.Meq+self.Menersys+self.propFFmass+self.Mvtolprop+self.Mfus+self.Mwing+self.Mhortail+self.Mverttail
        
    def getTailPars(self):
        self.xACh=self.xAC*self.MGC+self.xMGC+self.geodat.lh
        self.xLEMACh=self.xACh-0.25*self.MGCh
        self.xCrh=self.xLEMACh-np.tan(self.lambdaLEh)*self.yMGCh
        self.xcgh=self.xCrh+np.tan(self.lambdaLEh)*(0.38*self.bh/2)+0.42*self.getTailChord(0.38*self.bh/2)
        
    def getCircle(self, d, orig):
        self.loc=[[orig[0]+(d/2.)*np.cos(i) for i in np.linspace(0, 2*np.pi, 100)],
                         [orig[1]+(d/2.)*np.sin(i) for i in np.linspace(0, 2*np.pi, 100)]]
        return self.loc
        
    def getStability(self):
        self.getCoeffs()        
        self.xAC=self.acdat.Xac-(1.8/self.CLaAh)*((self.geodat.widthfus*self.geodat.hfus*self.xCr)/(self.S*self.MGC))+(0.273/(1+self.taper))*((self.geodat.widthfus*(self.S/self.b)*(self.b-self.geodat.widthfus))/((self.b*2.15*self.geodat.widthfus)*self.MGC**2))*np.tan(np.radians(self.acdat.lambdac4))
        
        
        self.one=(self.CLah/self.CLaAh)*(1-self.deda)*(self.geodat.lh/self.MGC)*self.VhV2
        self.two=(self.xAC-self.acdat.stabmarg)/self.one
        self.three=(self.CLh/self.CLAh)*(self.geodat.lh/self.MGC)*self.VhV2
        self.four=((self.CMac/self.CLAh)-self.xAC)/self.three
        self.xcgs=np.linspace(-2.1, 2.1, 1000)
        
        self.maneuverable=(1/self.one)*self.xcgs-self.two
        self.controllable=(1/self.three)*self.xcgs+self.four
        
    def getFFPropulsion(self):
        if self.xCr+self.rootchord>max(self.xcab):
            self.xffprop=self.xCr+self.rootchord
        else:
            self.xffprop=max(self.xcab)
    
    def getCoeffs(self):
        self.beta=np.sqrt(1-self.missiondat.Mach**2)
        self.lambdac2=np.arctan(0.5*self.b*np.tan(np.radians(self.lambdaLE))+(0.5*(self.tipchord-self.rootchord)/(0.5*self.b)))
        self.CLalphaw=(2*np.pi*self.A)/(2+np.sqrt(4+((self.A*self.beta/0.95)**2)*(1+(np.tan(self.lambdac2)**2)/(self.beta**2))))
        self.CLaAh=self.CLalphaw*(1+2.15*(self.geodat.widthfus/self.b))*(self.Swingw/self.S)+(np.pi/2)*((self.geodat.widthfus**2)/self.S)
    
        self.lambdac2h=np.arctan(0.5*self.bh*np.tan(np.radians(self.lambdaLEh))+(0.5*(self.tch-self.rch)/(0.5*self.bh)))
        self.CLah=(2*np.pi*self.taildat.Ah)/(2+np.sqrt(4+((self.taildat.Ah*self.beta/0.95)**2)*(1+(np.tan(self.lambdac2h)**2)/(self.beta**2))))
        
        self.rcoeff=self.geodat.lh/(0.5*self.b)
        self.mtv=0 #As we have tail-mounted propulsive system
        self.k1=((0.1124+0.1265*np.radians(self.acdat.lambdac4)+0.1766*np.radians(self.acdat.lambdac4)**2)/(self.rcoeff**2))+(0.1024/self.rcoeff)+2
        self.k2=0.1124/(self.rcoeff**2)+(0.1024/self.rcoeff)+2
        
        self.deda=(self.k1/self.k2)*((self.rcoeff**2/(self.rcoeff**2+self.mtv**2))*(0.4876/np.sqrt(self.rcoeff**2+0.6319+self.mtv**2))+(1+((self.rcoeff**2/(self.rcoeff**2+0.7915+5.0734*self.mtv**2))**0.3113))*(1-np.sqrt(self.mtv**2/(1+self.mtv**2))))*(self.CLalphaw/(np.pi*self.A))
        self.VhV2=self.acdat.VhV2
        
        self.CLh=-0.35*self.taildat.Ah**(1./3.)
        self.CLAh=self.acdat.CLAh
        self.CMacw=self.acdat.Cmac0*(self.A*(np.cos(np.radians(self.acdat.lambdac4))**2)/(self.A+2*np.cos(np.radians(self.acdat.lambdac4))))#Ask Melissa for an actual value of the airfoil Cmac 
        self.CMdeltaf=-1.8*(1-2.5*(self.geodat.widthfus/self.geodat.lfus))*((np.pi*self.geodat.widthfus*self.geodat.hfus*self.geodat.lfus)/(2*self.S*self.MGC))*(self.acdat.CL0/self.CLaAh)
        self.CMac=self.CMacw+self.CMdeltaf
    
    def calcFunctions(self):
        self.getPositions()
        self.addDiameters()
        self.getComponentWeights()
        self.getCabinGeometry()
        self.getWingGeometry()
        self.getFFPropulsion()
        self.getVTOLPropCG()
        self.getStability()
        self.getTailPars()
        self.getTailGeometry()
        self.getCGs()
       
    
    def drawPlanform(self):
        if self.show:
            plt.figure()
            plt.plot(self.xwing, self.ywing, label='Sweep={0} deg'.format(self.acdat.lambdac4), linewidth=0.7)
            plt.plot(self.xtail, self.ytail, label='Tail Sweep={:.02f} deg'.format(self.lambdac4h), linewidth=0.7)
            plt.plot([self.xMGC, self.xMGC+self.MGC], [self.yMGC, self.yMGC], label='Mean Geometric Chord', linewidth=0.5, linestyle='dashed')
            plt.plot(self.xcab, self.ycab, label='Cabin', color='k', linewidth=0.74)
            
            for i in range(self.Nprops):
                plt.plot(self.propdraw[i][0], self.propdraw[i][1], color='k', linewidth=0.74)
            if self.plotCGs:
                plt.scatter(self.xcgwing, 0, marker='o', s=20, label=r'$X_{cg_{wing}}$')
                plt.scatter(self.xcgfus, 0, marker='o', s=20, label=r'$X_{cg_{cabin}}$')
                plt.scatter(self.propxcg, 0, marker='o', s=20, label=r'$X_{cg_{prop}^{VTOL}}$')
                plt.scatter(self.empxcg, 0, marker='o', s=20, label=r'$X_{cg_{empennage}}$')
                plt.scatter(self.xffprop, 0, marker='o', s=20, label=r'$X_{cg_{prop}^{FF}}$')
                plt.scatter(self.xcgenersys, 0, marker='o', s=20, label=r'$X_{cg_{energy}}$')
                plt.scatter(self.xcgpl, 0, marker='o', s=20, label=r'$X_{cg_{PL}}$')
                plt.scatter(self.totalxcg, 0, marker='o', s=35, label=r'$X_{cg_{TOTAL}}$')
            plt.legend()
            plt.axes().set_aspect('equal', 'datalim')
            fig, ax1 = plt.subplots()
            ax1.plot(self.xcgs, self.controllable, label='Controllability curve', linewidth=4)
            ax1.plot(self.xcgs, self.maneuverable, label='Stability curve', linewidth=4)
            ax1.plot([self.mincg, self.maxcg], [self.ShS, self.ShS], label='CG range', linewidth=7)
            ax1.plot([ax1.get_xlim()[0], self.mincg], [self.ShS, self.ShS], 'k-', linewidth=0.5)
            plt.xlabel(r'$\frac{X_{cg}}{MAC}$, [-]', fontsize=18)
            plt.ylabel(r'$\frac{S_{h}}{S}$', fontsize=18)
            ax1.set_ylim(0)
#            ax1.set_xlim(-1.2, 1.2)
            ax1.fill_between(self.xcgs, y1=self.controllable, y2= 0, 
                 color="#fac2c2", zorder=0)
            ax1.fill_between(self.xcgs, y1=self.maneuverable, y2= 0, 
                 color="#fac2c2", zorder=0)   

            plt.legend()
            plt.show()
        
    def getLoadingDiag(self):
        plt.figure()
        sns.set(context='paper', font_scale=1.5, style='whitegrid')
        self.paxpos1=(self.geodat.lfus*np.array([self.geodat.Xcgpilot,
                                               self.geodat.Xcgseat3, self.geodat.Xcgpatient, 
                                               self.geodat.Xcgseat4])-self.xMGC)/self.MGC
        self.paxpos2=np.array(self.paxpos1[::-1])
        self.cgs1=np.array([(self.xcgoew-self.xMGC)/self.MGC])
        self.masses1=np.array([self.Moew])
        self.cgs2=np.array(self.cgs1)
        self.masses2=np.array(self.masses1)
        for j in range(len(self.paxpos1)):
            self.masses1=np.append(self.masses1, self.Moew+(j+1)*self.missiondat.Mpax)
            self.cgs1=np.append(self.cgs1, (self.cgs1[-1]*self.masses1[-2]+self.missiondat.Mpax*self.paxpos1[j])/(self.masses1[-1]))       
        for j in range(len(self.paxpos2)):
            self.masses2=np.append(self.masses2, self.Moew+(j+1)*self.missiondat.Mpax)
            self.cgs2=np.append(self.cgs2, (self.cgs2[-1]*self.masses2[-2]+self.missiondat.Mpax*self.paxpos2[j])/(self.masses2[-1]))       
        self.cgstotal=np.append(self.cgs1, self.cgs2)
        self.mincg=min(self.cgstotal)-0.05
        self.maxcg=max(self.cgstotal)+0.05
        self.tmp=np.intersect1d(self.xcgs[self.xcgs>=self.mincg], self.xcgs[self.xcgs<=self.maxcg])
        self.idx=[self.xcgs[i] in self.tmp for i in range(len(self.xcgs))]
        self.correct=np.append(self.controllable[self.idx], self.maneuverable[self.idx])
        self.ShS=max(self.correct[self.correct>0])
        if self.show:
            plt.plot(self.cgs1, self.masses1, '-D', label='PAX loading, front-to-back')
            plt.plot(self.cgs2, self.masses2, '-D', label='PAX loading, back-to-front')
            plt.xlabel(r'$\frac{X_{cg}}{MAC}$ [-]')
            plt.ylabel('M [kg]')
            plt.legend()
g=Geometry()
g.show=True
g.calcFunctions()
g.getLoadingDiag()
g.drawPlanform()
print('New Sh value: {} m2'.format(g.S*g.ShS))