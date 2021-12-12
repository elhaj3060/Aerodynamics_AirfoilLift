
import HSPM,geometryGenerator

# Creating the geometry and source panels
panels = geometryGenerator.GenerateNACA4digit(maxCamber=2.0,
                                              positionOfMaxCamber=4.0,
                                              thickness=12.0,
                                              pointsPerSurface=50)

# Instantiating HSPM class to compute the pressure solution on the given geometry
prob = HSPM.HSPM(listOfPanels = panels, alphaRange = [0.0,5.0,10.0,20.0,25.0])
# Solving...
prob.run()
# Extracting alpha max and cl max with Valarezo
alphaMax, clMax = prob.findAlphaMaxClMax(valarezoCriterion=14.0)
print('AlphaMax= %.2lf ClMax= %.4lf' % (alphaMax,clMax))