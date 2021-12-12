import numpy as np
from math import pi, sin, cos, sqrt, atan , log
import matplotlib.pyplot as plt

def cylinder_model(nb_panels, rayon):
    # Cette fonction permet de definir la geomtrie du cylindere

    # 
    mid_pts_angles = np.linspace(0,2*pi,nb_panels+1)
    mid_pts_angles = mid_pts_angles[:-1]
 
    edge_pts_angles = [(2*pi)/(nb_panels*2)+i for i in mid_pts_angles]
    

    edge_pts_x = [rayon*cos(i) for i in edge_pts_angles]
    edge_pts_y = [rayon*sin(i) for i in edge_pts_angles]

    # Coordonnee cartesienne des points du sommet
    edge_pts = zip(edge_pts_x,edge_pts_y)
    edge_pts = list(edge_pts)

    mid_pts_x = [(edge_pts_x[i]+edge_pts_x[i+1])/2 for i in range(0,len(edge_pts_x)-1)]
    mid_pts_y = [(edge_pts_y[i]+edge_pts_y[i+1])/2 for i in range(0,len(edge_pts_y)-1)]


    mid_pts_x.append((edge_pts_x[0]+edge_pts_x[-1])/2)
    mid_pts_y.append((edge_pts_y[0]+edge_pts_y[-1])/2)

    mid_pts_x = [round(i,6) for i in mid_pts_x]
    mid_pts_y = [round(i,6) for i in mid_pts_y]

    # Coordonnees cartesienne des points milieu

    mid_pts = zip(mid_pts_x,mid_pts_y)
    mid_pts = list(mid_pts)

    # mid_point[numero_du_point][x --> 0 ou y --> 1]

    
    edge_pts = np.array(edge_pts,dtype='float32')
    
    beta_angles = np.linspace(pi*2/nb_panels,2*pi,nb_panels,dtype='float32')
    
    phi_angles = []

    for i in beta_angles:
        if (i-pi/2) >= 0 :
            phi_angles.append(i-pi/2)
        else:
            phi_angles.append(i-pi/2 + 2*pi)



    return mid_pts,edge_pts,beta_angles,phi_angles,edge_pts_x,edge_pts_y


def integral_calculation(mid_pts,edge_pts,beta_angles,phi_angles):
    integrals_lambda = []
    integrals_vitesse = []
    
    '''
    Fonction pour calculer les integrales pour trouver les lambdas
    '''
    
    for i in range(0,len(mid_pts)):
        
        I_lambda = []
        I_vitesse = []
        
        for j in range(0,len(mid_pts)-1):

            A = (  -(mid_pts[i][0]-edge_pts[j+1][0])*cos(phi_angles[j])
                   -(mid_pts[i][1]-edge_pts[j+1][1])*sin(phi_angles[j])  )

            B = ((mid_pts[i][0] - edge_pts[j+1][0])**2 
                 + (mid_pts[i][1] - edge_pts[j+1][1])**2 )

            C = sin(phi_angles[i]-phi_angles[j])

            D = ( (mid_pts[i][1] - edge_pts[j+1][1])*cos(phi_angles[i])
                - (mid_pts[i][0] - edge_pts[j+1][0])*sin(phi_angles[i]) )

            S = ( sqrt(( edge_pts[j][0]- edge_pts[j+1][0])**2 
                      + (edge_pts[j][1] - edge_pts[j+1][1])**2) )

            E = sqrt(B-A**2)

            I_l = (C/2) * log((S**2 + 2*A*S+B)/B) + ((D-A*C)/E)*(atan((S+A)/E)-atan(A/E))
            
            I_v = ((D-A*C)/(2*E))*log((S**2+2*A*S+B)/B)-C*(atan((S+A)/E)-atan(A/E))

            I_lambda.append(I_l)
            I_vitesse.append(I_v)

        #La boucle est refaite pour connecter le dernier point et le premier pour completer la boucle 
        # Exemple : pts1 et pts2 / 2 et 3 /...../ 7et 8 ensuite il faut faire la derniere interation ----> 8 et 1


        A = (  -(mid_pts[i][0]-edge_pts[0][0])*cos(phi_angles[-1])
               -(mid_pts[i][1]-edge_pts[0][1])*sin(phi_angles[-1])  )

        B = (mid_pts[i][0] - edge_pts[0][0])**2 + (mid_pts[i][1] - edge_pts[0][1])**2 

        C = sin(phi_angles[i]-phi_angles[-1])

        D = ((mid_pts[i][1] - edge_pts[0][1])*cos(phi_angles[i])
            - (mid_pts[i][0] - edge_pts[0][0])*sin(phi_angles[i]))

        S = sqrt(( edge_pts[-1][0]- edge_pts[0][0])**2
            + (edge_pts[-1][1] - edge_pts[0][1])**2)

        E = sqrt(B-A**2)

        I_l = (C/2) * log((S**2 + 2*A*S+B)/B) + ((D-A*C)/E)*(atan((S+A)/E)-atan(A/E))
        
        I_v = ((D-A*C)/2*E)*log((S**2+2*A*S+B)/B)-C*(atan((S+A)/E)-atan(A/E))
        
        I_vitesse.append(I_v)

        I_lambda.append(I_l)

        integrals_lambda.append(I_lambda)
        integrals_vitesse.append(I_vitesse)


    return integrals_lambda, integrals_vitesse



def solveur_source_adim(integrals_lambda,beta_angles):
    '''
    Fonction pour calculer a partir des integrales et angles beta la grandeur
    des sources adimensionalisée
    '''


    for i in range(0,len(integrals_lambda)):
        for j in range(0,len(integrals_lambda)):
            if i == j:
                integrals_lambda[i][j] = pi

    A = np.array(integrals_lambda,dtype='float32')
    B = []
    for i in beta_angles:
        B.append(-cos(i))
    B = np.array(B,dtype='float32')

    
    adim_lambdas = np.linalg.solve(A,B)

    return adim_lambdas

def solveur_vitesses_adim(beta_angles, adim_lambdas,integrals_vitesse):
    
    adim_vitesse = []
    for i in range(len(integrals_vitesse[0])):
        somme_vitesse_adim = 0
        for j in range(len(integrals_vitesse[0])):

    	    somme_vitesse_adim += integrals_vitesse[i][j]*adim_lambdas[j]
        
        vitesse = sin(beta_angles[i]) + somme_vitesse_adim

        adim_vitesse.append(vitesse)


        
    return adim_vitesse
    


if __name__ == "__main__":

    mid_pts,edge_pts,beta_angles,phi_angles,edge_pts_x,edge_pts_y = cylinder_model(12,1)
    
    integrals_lambda,integrals_vitesse = integral_calculation(mid_pts,edge_pts,beta_angles,phi_angles)
    
    #--------------------------------------------------------------------------
    #                          Question 1a)
    #--------------------------------------------------------------------------
    
    adim_lambdas = solveur_source_adim(integrals_lambda,beta_angles)
    print('Le vecteur des sources adimensionaliée :', end='\n')
    print(np.round(adim_lambdas,3), end= '\n\n')
    
    #--------------------------------------------------------------------------
    #                           Question 1b)
    #--------------------------------------------------------------------------
    
    print('La somme des sources : ' ,end='')
    print(round(np.sum(adim_lambdas),2),end='\n')
    
    #--------------------------------------------------------------------------
    #                           Question 1c)
    #--------------------------------------------------------------------------
    
    adim_vitesse = solveur_vitesses_adim(beta_angles, adim_lambdas, integrals_vitesse)
    cp = []
    for i in range(0,len(adim_vitesse)):
        cp.append(1 - adim_vitesse[i]**2)
        
    plt.plot(beta_angles,cp,'x',beta_angles,[1-4*sin(i)**2 for i in beta_angles])
    plt.xlabel("Angle beta (rads)")
    plt.ylabel("Cp")
    plt.legend(('sol_numerique','sol_theorique'))
    plt.show()
    
    #--------------------------------------------------------------------------
    #                           Question 1d)
    #--------------------------------------------------------------------------
    
    edge_pts_x.append(edge_pts_x[0])
    edge_pts_y.append(edge_pts_y[0])
    panel_length = []
    for i in range(len(edge_pts_x)-1):
        panel_length.append( sqrt( ((edge_pts_x[i+1]-edge_pts_x[i])**2)+((edge_pts_y[i+1]-edge_pts_y[i])**2) ))
        
        
    forces = np.empty((len(panel_length)),dtype='float32')
    print("\nForces : \n")
    for i in range(len(panel_length)):
        forces[i] = (cp[i]*panel_length[i])
        print("Pour beta de " + str(round(np.degrees(beta_angles[i]),2)) + " deg : " + str(forces[i])) 
        
        
    
        
    
