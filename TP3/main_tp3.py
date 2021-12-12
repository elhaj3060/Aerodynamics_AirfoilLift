import sys

sys.path.append("./VLM_etudiants")
from VLM_etudiants import vlm
import matplotlib.pyplot as plt
import numpy as np
import os
from math import exp, pi
from LP import Methode_LP
## Question1 a)
span_inf = 10000
prob1 = vlm.VLM(ni=5,
                nj=50,
                chordRoot=1.0,
                chordTip=1.0,
                twistRoot=0.0,
                twistTip=0.0,
                span=span_inf,
                sweep=0.0,
                Sref=span_inf,
                referencePoint=[0, 0, 0],
                wingType=1,
                alphaRange=[2, 4],
                )
prob1.run()
X = np.array(prob1.alphaRange) * 0.01745
Y = np.array(prob1.CL)
print(np.polyfit(X, Y, 1)[0]);
# Resultat : 6.215053375802405

## Question 1b)

def sref_b(span,chord_root,chord_tip):
    return (chord_root + chord_tip)*span/2


sweep = [0,30,45,60]
spans = [0.5,1,2,4]
Cl_alpha_list_b = []
AR = [2*i**2/sref_b(i,1,1) for i in spans]

for j in sweep:
    for i in spans:
        prob2 = vlm.VLM(ni=5,
                        nj=40,
                        chordRoot=1.0,
                        chordTip=1.0,
                        twistRoot=0.0,
                        twistTip=0.0,
                        span=i,
                        sweep=j,
                        Sref=sref_b(i,1,1),
                        referencePoint=[0, 0, 0],
                        wingType=1,
                        alphaRange=[0,5],
                        )

        prob2.run()
        # Calcul de la pente cl_alpha
        Cl_alpha_list_b.append((prob2.CL[1]-prob2.CL[0])/((prob2.alphaRange[1]-prob2.alphaRange[0])*0.01745))

    plt.plot(AR, Cl_alpha_list_b)
    Cl_alpha_list_b.clear()
plt.legend(("\u039B : 0 deg"," \u039B :30 deg", "\u039B : 45 deg","\u039B : 60 deg"))
plt.xlabel("AR")
plt.ylabel("CL_")
plt.show()

## Question 1c) -- Generation des fichiers

def span_c(AR,chord_root,chord_tip):
    return AR*(chord_root + chord_tip)/4
def sref_c(span,chord_root,chord_tip):
    return (chord_root + chord_tip)*span/2

AR = 4 # Dans le manuel
Cl_list_c = []
alpha_c = 4
sweep = [0,45,135]
for i in sweep:
    prob3 = vlm.VLM(ni=5,
                    nj=50,
                    chordRoot=1.0,
                    chordTip=1.0,
                    twistRoot=0.0,
                    twistTip=0.0,
                    span=span_c(AR,1,1),
                    sweep=i,
                    Sref=sref_c(span_c(AR,1,1),1,1),
                    referencePoint=[0, 0, 0],
                    wingType=1,
                    alphaRange=[alpha_c],
                    )


    prob3.run()
    donnees_cl = np.loadtxt('Spanload_A' + str(alpha_c) + '.00.dat', skiprows=1)
    corde,dist_cl = donnees_cl[:, 0],donnees_cl[:, 1]
    ratio_cl = dist_cl / prob3.CL
    plt.plot(corde / span_c(AR,1,1), ratio_cl)
plt.ylabel("cl/CL")
plt.xlabel("% Demi-envergure")
plt.legend(("\u039B : 0 deg"," \u039B :45 deg", "\u039B : 135 deg"))
plt.show()

## Question 1d) -- Generation des fichiers
def span_d(AR,chord_root,chord_tip):
    return AR*(chord_root + chord_tip)/4
def sref_d(span,chord_root,chord_tip):
    return (chord_root + chord_tip)*span/2

alpha_d = 4
Cl_list_d = []
AR = 7.28
tip = [0,0.4,0.6,1]
for i in tip:
    prob4 = vlm.VLM(ni=5,
                    nj=40,
                    chordRoot=1.0,
                    chordTip=i,
                    twistRoot=0.0,
                    twistTip=0.0,
                    span=span_d(AR,1,i),
                    sweep=0,
                    Sref=sref_d(span_d(AR,1,i),1,i),
                    referencePoint=[0, 0, 0],
                    wingType=1,
                    alphaRange=[alpha_d],
                    )

    prob4.run()

    donnees_cl = np.loadtxt('Spanload_A' + str(alpha_d) + '.00.dat', skiprows=1)
    corde, dist_cl = donnees_cl[:, 0], donnees_cl[:, 1]
    ratio_cl = dist_cl / prob4.CL
    ratio_corde = corde / span_d(AR, 1, i)
    plt.plot(ratio_corde, ratio_cl)
plt.ylabel("cl/CL")
plt.ylim(0,1.4)
plt.xlabel("% Demi-envergure")
plt.legend(("\u03BB : 0 ", " \u03BB :0.4 ", "\u03BB : 0.6 ","\u03BB : 1"))
plt.show()
## Question 1e) -- Generation des fichiers
span_e = 10
alpha_e = [0,5]
sref_e =  0.5*pi*span_e*0.5 # chord_root = chord_tip =  1
AR = ((span_e*2)**2)/(2*sref_e)
prob5 = vlm.VLM(ni=5,
                nj=40,
                chordRoot=1.0,
                chordTip=1,
                twistRoot=0.0,
                twistTip=0.0,
                span=span_e,
                sweep=0,
                Sref=sref_e,
                referencePoint=[0, 0, 0],
                wingType=3,
                alphaRange=alpha_e,
                )


prob5.run()
# Calcul de la pente cl_alpha
pente_cl_alpha_e = ((prob5.CL[1]-prob5.CL[0])/((prob5.alphaRange[1]-prob5.alphaRange[0])*0.01745))
# calcul de la pente
donnees_cl = np.loadtxt('Spanload_A' + str(alpha_e[1]) + '.00.dat', skiprows=1)
corde, dist_cl = donnees_cl[:, 0], donnees_cl[:, 1]
c = 1*(1-(corde/span_e)**2)**0.5
c_avg = sref_e/span_e


plt.plot(corde/span_e, (dist_cl/prob5.CL[1])*(c/c_avg))
plt.ylabel("CL")
plt.xlabel("Demi-envergure")
plt.show()

oswald = prob5.CL[1]**2/(pi*AR*prob5.CD[1])
print("coefficient de Oswald pour alpha de " + str(alpha_e[1]) +  " :  " + str(oswald))

pente_enonce = (2*pi)/(1+2/(AR))
print("pente theorique : " + str(pente_enonce))
print("pente calculee : " + str(pente_cl_alpha_e))


## Question 2
alpha_2 = [5]
sref_e =  0.5*pi*0.5 # chord_root = chord_tip =  1
AR = ((1*2)**2)/(sref_e*2)
prob_LP = vlm.VLM(ni=5,
                nj=50,
                chordRoot=1.0,
                chordTip=1,
                twistRoot=0.0,
                twistTip=0.0,
                span=1,
                sweep=0,
                Sref=sref_e,
                referencePoint=[0, 0, 0],
                wingType=3,
                alphaRange=alpha_2,
                )


prob_LP.run()
#pente_cl_alpha_e = ((prob5.CL[1]-prob5.CL[0])/((prob5.alphaRange[1]-prob5.alphaRange[0])*0.01745))
# calcul de la pente
donnees_cl = np.loadtxt('Spanload_A' + str(alpha_2[0 ]) + '.00.dat', skiprows=1)
corde, dist_cl = donnees_cl[:, 0], donnees_cl[:, 1]




CL_LP, CDi_LP, oswald_LP = Methode_LP(corde,dist_cl,1,AR,1,50)
print("Oswald : " + str(oswald_LP))
print("CL_LP : " + str(CL_LP))
print("CDi_LP : " + str(CDi_LP))

    ## Question 3
    # Apres plusieurs essais
    # twist -4 > angle 7.8
    # twist 0 > angle 6.2
    # twist 4 > angle 4.6
    tip_twist = [-4,0,4]
    span_3 = 5
    chord_tip3 = 0.4
    sref_3 =  (1+chord_tip3)*span_3*0.5
    AR = ((span_3*2)**2)/(sref_3*2)
    alpha_3 = [[7.8],[6.2],[4.6]]
    for twist,alpha in zip(tip_twist,alpha_3):
        prob6 = vlm.VLM(ni=5,
                        nj=50,
                        chordRoot=1.0,
                        chordTip=chord_tip3,
                        twistRoot=0.0,
                        twistTip=twist,
                        span=span_3,
                        sweep=35,
                        Sref=sref_3,
                        referencePoint=[0, 0, 0],
                        wingType=2,
                        alphaRange=alpha,
                        )
        prob6.run()
        donnees_cl = np.loadtxt('Spanload_A' + str(alpha[0]) + '0.dat', skiprows=1)
        list_corde, dist_cl = donnees_cl[:, 0], donnees_cl[:, 1]

        c = 1*(1-(list_corde/span_3)**2)**0.5
        c_avg = sref_3/span_3
        ratio_cl = [dist_cl[i]*c[i]/c_avg for i in range(len(c))]


        plt.plot(list_corde/span_3, ratio_cl)
        CL_LP_3, CDi_LP_3, oswald_LP_3 = Methode_LP(list_corde, dist_cl, 1, AR, 5, 50)
        print("Oswald pour 3 : " + str(oswald_LP_3))
        print("CDi_LP pour 3: " + str(CDi_LP_3))
        print(CL_LP_3)



    plt.ylabel("cl/CL")
    plt.xlabel("% Demi-envergure")
    plt.legend(("twist: -4 ", " twist : 0 ", "twist : 4 "))
    plt.show()

