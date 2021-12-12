# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 23:22:56 2020

@author: rouga
"""
from libblasius import integrateTrapeze
def solveur_thwaites_turbulent(ue,x,uinf,i_pour_transition):
    pos = 0
    while  pos <= i_pour_transition:
        theta = (integrateTrapeze(x,ue**5,0,pos)*0.45)/(ue[pos]**6)
        lmbda =(theta)*(ue[pos+1]-ue[pos])/(x[pos+1]-x[pos])
        pos += 1
    
    cf = True
    H=1.4
    
    i_pour_transition += 1
    while H < 3.0 and cf>0.00 and i_pour_transition <= len(ue)-pos-1: 
    
        
        cf = 0.246*(10**(-0.678*H))*(ue[i_pour_transition]*theta)**(-0.268)
        
        dx = x[i_pour_transition+1]-x[i_pour_transition]
        
        d_ue_dx = ( ue[i_pour_transition+1]-ue[i_pour_transition])/(dx)
        
        
        theta = ((cf/2) - (H+2)*(theta/ue[i_pour_transition])*d_ue_dx)*dx + theta
                
        H = (((-H*(H-1)*(3*H-1)*(theta/ue[i_pour_transition])*d_ue_dx + H*(3*H-1)*(cf/2)-((3*H-1)**2)*(0.0056/2)*(ue[i_pour_transition]*theta)**(-1/6))/theta)*dx)+H
        i_pour_transition += 1
        
    print("Theta : "+ str(theta) )
    print("cf : " + str(cf))
    print("H :" + str(H))    
    return i_pour_transition
    
    