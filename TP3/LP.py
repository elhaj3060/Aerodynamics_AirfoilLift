import numpy as np
from math import *

def Methode_LP(liste_Coord,liste_cl_locale,chord_root,AR,semi_span,nb_panels_y):

    liste_Coord = np.array(liste_Coord)
    liste_cl_locale= np.array(liste_cl_locale)



    theta = np.arccos(-1*liste_Coord/ semi_span)

    cy = (((chord_root)**2)*(1-(liste_Coord/semi_span)**2))**0.5


    gamma =  0.5*np.multiply(liste_cl_locale,cy)

    nbrs_imp = [i for i in range(1,nb_panels_y*2,2)]
    
    
    A = np.zeros((len(theta),len(nbrs_imp)))
    for i in range(len(theta)):
        for j in range(len(nbrs_imp)):
            A[i,j] = sin(theta[i]*nbrs_imp[j])*4*1


    coeffs_A = np.linalg.solve(A,gamma)

    Del = 0
    
    
    for i in range(1,nb_panels_y):
        Del_tmp = ( nbrs_imp[i] * (coeffs_A[i]/coeffs_A[0])**2 )
        Del += Del_tmp

    oswald = 1/(1+Del)
    
    CL = pi*coeffs_A[0]*AR
    
    CDi = CL**2/(pi*AR*oswald)
    
    return CL, CDi , oswald



