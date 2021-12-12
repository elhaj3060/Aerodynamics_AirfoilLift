import sys
import matplotlib.pyplot as plt
sys.path.append("./HSPM_etudiants")
import  HSPM
import geometryGenerator as GeoGen


panels = GeoGen.ReadPoints("CW_NACA1410.dat")
GeoGen.ReadPoints("CW_NACA1410.dat")
alpha_att = [i for i in range(-20,21,5)]

if __name__ == "__main__" :

    prob= HSPM.HSPM(listOfPanels = panels,alphaRange = alpha_att)


    prob.run()
   
    # Question2 a- 
    plt.figure()
    
    plt.subplot(211)
    plt.plot(alpha_att,prob.CL)
    plt.title("Alpha vs CL (NACA1410)")
    
    plt.subplot(212)
    plt.plot(alpha_att,prob.CM)
    plt.title("Alpha vs CM (NACA1410)")
    plt.show()