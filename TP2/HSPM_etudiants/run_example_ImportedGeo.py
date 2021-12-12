import HSPM,geometryGenerator

# Creating the geometry and source panels
panels = geometryGenerator.ReadPoints("airfoil_points.dat")

# Instantiating HSPM class to compute the pressure solution on the given geometry
prob = HSPM.HSPM(listOfPanels = panels, alphaRange = [0.0,5.0,10.0,20.0,25.0])
# Solving...
prob.run()
# Extracting alpha max and cl max with Valarezo
alphaMax, clMax = prob.findAlphaMaxClMax(valarezoCriterion=14.0)
print('AlphaMax= %.2lf ClMax= %.4lf' % (alphaMax,clMax))