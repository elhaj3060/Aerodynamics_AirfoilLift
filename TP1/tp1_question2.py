import sys
import matplotlib.pyplot as plt
sys.path.append("./HSPM_etudiants")
import  HSPM
import geometryGenerator as GeoGen


panels = GeoGen.ReadPoints("CW_NACA6412.dat")
alpha_att = [i for i in range(-20,21,5)]
prob_a = HSPM.HSPM(listOfPanels = panels,alphaRange = alpha_att)
prob_c = HSPM.HSPM(listOfPanels = panels,alphaRange = alpha_att,
                   referencePoint=[0.25,0,0])

if __name__ == "__main__" :
    prob_a.run()
   
    # Question2 a- 
    plt.figure()
    
    plt.subplot(311)
    plt.plot(alpha_att,prob_a.CL)
    plt.title("Alpha vs CL (NACA6412)")
    
    plt.subplot(312)
    plt.plot(alpha_att,prob_a.CD)
    plt.title("Alpha vs CD (NACA6412)")
    
    plt.subplot(313)
    plt.plot(alpha_att,prob_a.CM)
    plt.title("Alpha vs CM (NACA6412)")
    plt.show()
    # Question2 c- 
    prob_c.run()
    plt.plot(alpha_att,prob_c.CM)
    plt.title("Alpha vs CM avec le point de reference au quart de la corde")
    plt.show()
    
    