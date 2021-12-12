import sys
import numpy as np
sys.path.append("./TD2_donnees")
sys.path.append("./HSPM_etudiants")

from libblasius import Blasius, integrateTrapeze
from SOLVEUR_THWAITES import solveur_thwaites
import matplotlib.pyplot as plt
import HSPM
from Vector3 import Vector3
from sourcePanel import sourcePanel
from tp1_question1 import cylinder_model
import math
from scipy.spatial import distance
from SOLVEUR_THWAITES_TURBULENT import solveur_thwaites_turbulent




# %% Question 1 - b i
def run_q1():
    print("\n\n --Question 1-- \n\n")
    couche_blasius = Blasius(N=101, numIer=20, x=1.0,
                             etaStart=0, etaEnd=10.0,
                             VInf = 1,tol=-16.0)
    couche_blasius.run()
    eta_h_cl = np.interp(0.995,couche_blasius.fp,couche_blasius.eta)
    x = np.linspace(0,1,101)
    xx = np.linspace(0,10,101)
    # %% Question 2 - b ii
    del_star = []
    for i in range(0,len(x)):
        del_star_i =integrateTrapeze(xx,1-couche_blasius.fp,0,i)
        del_star.append(del_star_i)
    
    # %% Question2 - b -iii
    theta = []
    for i in range(0,len(x)):
        theta_i =integrateTrapeze(xx,(1-couche_blasius.fp)*couche_blasius.fp,0,i)
        theta.append(theta_i)
    
    
# %% Question 2
        
def run_q2():
    print("\n\n --Question 2-- \n\n")
    
    uinf = 1
    c= 1                                                                                                                                                                                
    nb_pts = 4000
    lmbda_sep = -0.09
    x = np.linspace(0,c,nb_pts+1)
    ue = uinf*(1-x/c)
    
    sep,pos,lmbda_list = solveur_thwaites(ue,x,uinf,-0.09)
    print(sep,pos)   
    
    plt.plot([x[i]/c for i in range(0,pos-1)],lmbda_list)
    plt.ylabel("lambda")
    plt.xlabel("x\c")
    plt.title("Lambda en fonction de la position le long de la corde")
    plt.show()
 

# %% Donnees pour Q3 et Q4
nbr_panneaux = 150
data_tmp = cylinder_model(nbr_panneaux,1)[1]

# Coordonnees du profil
data = np.insert(data_tmp,0,data_tmp[-1],axis=0)
data = np.flipud(data)

# transformer les coordonnees en input lisible par HSPM (copier-coller de la fontion
# readPoints de geometryGenerator)
panels = []
for i in range(0,len(data)-1):
    pt     = data[i,:]
    ptNext = data[i+1,:]
    p1 = Vector3(pt[0],     0.0, pt[1])
    p2 = Vector3(ptNext[0], 0.0, ptNext[1])
    panels.append(sourcePanel(p1,p2))
    
    
    
model = HSPM.HSPM(listOfPanels=panels,alphaRange=[0])
model.run()

extrados_V = model.getUpperVtangential()[1]

extrados_V = np.array(extrados_V)
angles = np.linspace(0,180,len(extrados_V))

angles_rad = np.linspace(0,math.pi,len(extrados_V))   

s = [0]
dist_total = 0
for i in range(0,len(extrados_V)-1):
    dist = distance.cdist([data[i]],[data[i+1]])[0][0]
    dist_total += dist
    s.append(dist_total) 
    
# %% Question 3
    
def run_q3():
    print("\n\n --Question 3-- \n\n")
        
    
    point_sep,pos = solveur_thwaites(extrados_V,s,1,-0.09)[0:2]
    print("Separation laminaire : " + str(angles[pos]))
    
# %% Question 4
    
def run_q4():
    
    print("\n\n --Question 4-- \n\n")
    
    i_pour_10deg = 0
    while angles[i_pour_10deg] <= 10:
        i_pour_10deg +=1

    pos_turbulent = solveur_thwaites_turbulent(extrados_V,s,1,i_pour_10deg+1)
    print("Separation turbulent : " + str(angles[pos_turbulent]))

    
# %% les fonctions dans le main  
if __name__ == "__main__":

    run_q4()
    
    