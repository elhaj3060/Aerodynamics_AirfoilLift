# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 18:55:20 2020
"""

import sys
sys.path.append("./TD2_donnees")
import numpy as np

from libblasius import integrateTrapeze

def solveur_thwaites(ue,x,uinf,lmbda_crit):
    
    pos = 1
    lmbda = 5
    lmbda_list = []
    while lmbda > lmbda_crit and pos < len(x)-1:

        theta = (integrateTrapeze(x,ue**5,0,pos)*0.45)/(ue[pos]**6)
        lmbda =(theta)*(ue[pos+1]-ue[pos])/(x[pos+1]-x[pos])
        lmbda_list.append(lmbda)
        pos += 1

    return x[pos],pos, lmbda_list


    
   
    




