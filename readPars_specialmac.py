# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:29:56 2019

@author: chris (Miks edit)
"""

import os, pandas
from argparse import Namespace 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def load_input(Vars=True, sheet=None, namespacename='inp',name=None, info=None):
    if name == None:
        name = 'Parameters'
    if sheet == None:
        raise Exception('Specify the sheet name!')
    else:
        path = os.getcwd()
        data = pandas.read_excel(io=path+'/'+name+'.xlsx', sheet_name=sheet)
        dictvalues = dict(zip(data.Parameter, data.Value))
        if info==None:
            if Vars==True:
                inp = Namespace(**dictvalues)
                return inp
            else:
                return dictvalues
        else:
            return data

def plot_graph(x, y, xlabel=None, ylabel=None, newfig=True, label=None):
    sns.set(context='paper', font_scale=1.5, style='whitegrid')
    if type(x)==list or type(x)==np.matrixlib.defmatrix.matrix:
        x=np.array(x)
    if type(y)==list or type(y)==np.matrixlib.defmatrix.matrix:
        y=np.array(y)
    if xlabel==None:
        raise Exception('Indicate label for the x axis!')
    if ylabel==None:
        raise Exception('Indicate label for the y axis!')
    if newfig:
        plt.figure()
    if len(np.shape(x))==2:
        for i in range(np.shape(x)[0]):
            if label!=None:
                plt.plot(x[i], y[i], label=label[i])
            else:
                plt.plot(x[i], y[i])
    else:
        if label!=None:
                plt.plot(x, y, label=label)
        else:
                plt.plot(x, y)
        
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()#, loc='upper left')
    plt.show()
    
""""
Example use:

In []: inputdata = load_input()
In []: inputdata.e
Out[]: 0.8

"""