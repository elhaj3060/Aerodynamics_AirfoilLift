"******************************************************************************************************"
""" ***** TD2 - AER8270 : AÃ©rodynamique, Feng Yang Chen (1893644) & Baptiste Langlet (1840882) ***** """
"******************************************************************************************************"

#%% QUESTION 2
import numpy as np
from libblasius import integrateTrapeze

U_inf=100
nu = 1
s = np.linspace(0,1,1000)
U_e=U_inf*(1-s)
list(s)
list(U_e)
S=np.array(s)
U_E=np.array(U_e)


def ThwaitesNum(S, U_E):
    l=0
    i=0
    iteration = 1
    while (l >= -0.09) | (iteration < 100):
        integral = 0.45*integrateTrapeze(S,U_E**5,0,i)
        l = -U_inf*integral/(U_E[i+1]**6)
        i = i+1
        iteration = iteration + 1
    return l, i

l, iteration = ThwaitesNum(S, U_E)

#%% QUESTION 3
import sys
sys.path.append("./HSPM_etudiants")
import numpy as np
import HSPM
import math


# Du TD1 pour les coordonees d'un cylindre
def intI(n):
    X = []
    Y = []
    x = []
    y = []
    Phi = []
    
    
    # Longueur des panneaux avec un cercle unitaire
    S = math.sqrt(2*(1-math.cos(math.radians(360/n))))
    
    
    # Pour le premier angle différent
    phi_1 = math.radians(90)
    beta_1 = math.radians(-360/(2*n))
    X_1 = -math.cos(beta_1)
    Y_1 = math.sin(beta_1)
    x_1 = X_1 + (S/2)*math.cos(phi_1)
    y_1 = Y_1 + (S/2)*math.sin(phi_1)
    
    # Rajouter les premiers points à la liste de réponse
    X.append(X_1)
    Y.append(Y_1)
    x.append(x_1)
    y.append(y_1)
    Phi.append(phi_1)
    
    
    # Répéter les mêmes étapes pour les autres points 
    for i in range(2,n+1):
        phi = math.radians(90 - (360/n)*(i-1))
        beta = math.radians(-360/(2*n) + (360/n)*(i-1))
        X_i = -math.cos(beta)
        Y_i = math.sin(beta)
        x_i = X_i + (S/2)*math.cos(phi)
        y_i = Y_i + (S/2)*math.sin(phi)
        
        # Rajouter les résultats aux listes respectives
        X.append(X_i)
        Y.append(Y_i)
        x.append(x_i)
        y.append(y_i)
        Phi.append(phi)
    Theta = []
    for i in range(1,n+1):   
        Theta.append(math.radians(180-(360/n)*(i-1))+(np.pi))
        
    return X, Y, x, y, Phi, Theta 

a = intI(1000)

n=1000
# Write in txt file
file = open("HSPM_input.txt","a") 
for i in range(0,n):
    file.writelines([str(a[2][i-int(n/2)]) + ' ' + str(a[3][i-int(n/2)]) + '\n'])
file.writelines([str(a[2][-int(n/2)]) + ' ' + str(a[3][-int(n/2)]) ])
file.close()

#%% Lire le fichier contenant les coordonnees du cylindre
import geometryGenerator as geoGen

panels = geoGen.ReadPoints('HSPM_input.txt')
# Appeler HSPM pour trouver les vitesses
prob = HSPM.HSPM(listOfPanels = panels, alphaRange = [0])
prob.run()

# Obtenir les vitesses du HSPM
(lowerCoor, lowerV) = prob.getLowerVtangential()
(upperCoor, upperV) = prob.getUpperVtangential()

#%% Associer le bon theta aux bonnes vitesses

angle = []
for i in range(int(n/2)):
    angle.append(2*np.pi*i/n)
angle.append(np.pi)

# Thwaites
nu = 1
i=1
step = angle[1]-angle[0]
Theta = np.zeros((1,len(upperV)))
dU_e = np.zeros((1,len(upperV)))
lmba = np.zeros((1, len(upperV)))

# CALCUL DE LA DERIVEE
for i in range(len(upperV)-1):
		if i==0:
			dU_e[0,i]=(upperV[i+1]-upperV[i])/step
		elif i==(len(upperV)-1):
			dU_e[0,i]=(upperV[i]-upperV[i-1])/step
		else:
			dU_e[0,i]=(upperV[i+1]-upperV[i-1])/(2*step)

# Calcul pour les theta
for j in range(len(upperV)-2):
	Theta[0,j+1]=math.sqrt((nu/(upperV[j+1]**6))*(0.45*((step)/2)*(upperV[j]**5 + upperV[j+1]**5) + (((Theta[0,j]**2)*upperV[j]**6)/nu)))
# Calcul pour les lambdas
for k in range(len(upperV)-1):
    lmba[0,k]=((Theta[0,k]**2)/nu)*dU_e[0,k]
   
b = 1
item = 0    
while b :
    if lmba[0,item] > -0.09:
        item = item + 1
    else:
        b = 0
print(lmba[0,item])
print(item)
print(angle[item])

#%% QUESTION 4

nu = 1.568*(10**(-5))
# Separation entre laminaire et turbulent
lim_lam = round(10/180*len(upperV))

#Nombre de points turbulents
nbre_pnt_turb = len(upperV) - lim_lam

# Initialisation des parametres turbulents
Theta_turb = np.zeros((1,nbre_pnt_turb))
H = np.zeros((1,nbre_pnt_turb))
c_f = np.zeros((1,nbre_pnt_turb))

# Imposition des valeurs initiales
Theta_turb[0,0] = Theta[0,lim_lam]
H[0,0] = 1.4
c_f[0,0] = 0.246*10**(-0.678*H[0,0])*(Theta_turb[0,0]*upperV[lim_lam]/nu)**(-0.268)

# Refaire pour le reste
for i in range(nbre_pnt_turb-1):
    c_f[0,i+1] = 0.246*(10**(-0.678*H[0,i]))*((Theta_turb[0,i]*upperV[lim_lam+i]/nu)**(-0.268))
    dtheta = c_f[0,i]/2 - (H[0,i]+2)*Theta_turb[0,i]/upperV[lim_lam+i]*(dU_e[0,lim_lam+i])
    Theta_turb[0,i+1] = dtheta*step+Theta_turb[0,i]
    #
    dH = (-H[0,i]*(H[0,i]-1)*(3*H[0,i]-1)*Theta_turb[0,i]/upperV[lim_lam+i]*(dU_e[0,lim_lam+i]) + H[0,i]*(3*H[0,i]-1)*c_f[0,i]/2 - ((3*H[0,i]-1)**2)*0.0056*0.5*((Theta_turb[0,i]*upperV[lim_lam+i]/nu)**(-1/6)))/(Theta_turb[0,i+1])
    H[0,i+1] = dH*step+H[0,i]
    #
    # c_f[0,i+1] = 0.246*(10**(-0.678*H[0,i]))*((Theta_turb[0,i]*upperV[lim_lam+i]/nu)**(-0.268))
    # cf_turbulent[0,i+1]=0.246*10**(-0.678*H_turbulent[0,i])*(Ue[i+indice_turbulent]*theta_turbulent[0,i]/nu)**(-0.268)
		
# Condition d'arret
b = 1
item = 0    
while b :
    if (c_f[0,item] <= 0) | (H[0,item] >= 3) :
        b=0
    else:
        item = item + 1

# Resultat
print(c_f[0,item])
print(H[0,item])
print(angle[item+lim_lam])
