
import HSPM,geometryGenerator
import numpy as np

# Creating the geometry and source panels
panels = geometryGenerator.GenerateNACA4digit(maxCamber=0.0,
                                              positionOfMaxCamber=0.0,
                                              thickness=12.0,
                                              pointsPerSurface=120)

# Instantiating HSPM class to compute the pressure solution on the given geometry
prob = HSPM.HSPM(listOfPanels = panels, alphaRange = [0.0, 5.0, 10, 15])
# Solving...
prob.run()


print(prob.deltaCPvalarezo)

(alphaMax, clMax)=prob.findAlphaMaxClMax(valarezoCriterion=5.0)
print("alphaMax: %.3lf, clMax: %.3lf" % (alphaMax, clMax))

