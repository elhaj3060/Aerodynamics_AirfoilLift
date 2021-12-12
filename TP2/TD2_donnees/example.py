from libblasius import *

test = Blasius(	N=20, 
				numIter=10, 
				x=1.0,
				etaStart=0.0,
				etaEnd=5.0,
				VInf=5.0,
				tol=-16.0)

test.initLinear()
test.solveBlasius()
test.printResults()