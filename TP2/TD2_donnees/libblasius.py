# Standard libraries
import numpy as np
import math

# Methods and class for solving Blasius equation on flat plate

# Trapeze integration
def integrateTrapeze(xArray, yArray, iStart, iEnd):
   N  = min(len(xArray),len(yArray))

   if (iStart == iEnd) | (N == 1):
      return 0.0
   elif (iStart > N) | (iStart < 0) | (iEnd > N) | (iEnd < 0):
      print("Wrong starting or ending index for integration!")
      exit()
   else:
      result   = 0.0
      factor   = 1.0
      if (iStart > iEnd):
         factor   = -1.0
         tmp      = iEnd
         iEnd     = iStart
         iStart   = tmp


      for i in range(iStart+1, iEnd+1):
         result   += 0.5 * (yArray[i] + yArray[i - 1]) * (xArray[i] - xArray[i - 1])

      return (factor * result)

# Blasius class
class Blasius:
   # Initialization of the class
   def __init__(self, **kwargs):

      self.numNodes     = kwargs.pop('N',100)
      self.numIter      = kwargs.pop('numIter',100)
      self.x            = kwargs.pop('x',1.0)
      self.etaStart     = kwargs.pop('etaStart',0.0)
      self.etaEnd       = kwargs.pop('etaEnd',10.0)
      self.VInf         = kwargs.pop('VInf',1.0)
      self.tol          = kwargs.pop('tol',-16.0)
      self.rms          = 0.0
      self.rms0         = 0.0
      
      # Integration variables
      N                 = self.numNodes
      self.eta          = np.zeros(N)
      self.f            = np.zeros(N)
      self.fp           = np.zeros(N)
      self.fp0          = np.zeros(N)
      self.fpp          = np.zeros(N)
      self.F            = np.zeros(N)

      # Discretization 
      deta              = (self.etaEnd - self.etaStart) / (float(N - 1))

      for i in range(0,N):
         self.eta[i]    = self.etaStart + float(i) * deta

   # Initialization of the solution to a linear approximation
   def initLinear(self):

      N                 = self.numNodes

      for i in range(0,N):
         self.F[i]      = 0.0
         self.f[i]      = 0.0
         self.fp[i]     = self.eta[i] / self.etaEnd
         self.fp0[i]    = self.fp[i]
         self.fpp[i]    = 0.0

      # Boundary conditions
      self.f[0]         = 0.0
      self.fp[0]        = 0.0
      self.fp[N - 1]    = 1.0

   # Solving Blasius equation
   def solveBlasius(self):

      iterMax  = self.numIter
      N        = self.numNodes

      for it in range(0,iterMax):
         # Step 1: Saving fp in fp0
         for i in range(0,N):
            self.fp0[i]    = self.fp[i]

         # Step 2: Integration of fp to obtain f
         for i in range(0,N):
            self.f[i]      = integrateTrapeze(self.eta, self.fp, 0, i)

         # Step 3: Imposing boundary condition for f
         self.f[0]         = 0.0

         # Step 4: Integration of f to obtain F
         for i in range(0,N):
            self.F[i]      = integrateTrapeze(self.eta, self.f, 0, i)

         # Step 5: Imposing boundary condition for fpp
         self.fpp[0]       = (self.fp[1] - self.fp[0]) / (self.eta[1] - self.eta[0])

         # Step 6: Computing fpp from analytical relation
         for i in range(0,N):
            self.fpp[i]    = self.fpp[0] * math.exp(-0.5 * self.F[i])

         # Step 7: Integration of fpp to obtain fp
         for i in range(0,N):
            self.fp[i]     = integrateTrapeze(self.eta, self.fpp, 0, i)

         # Step 8: Imposing boundary condition for fp and computing rms for convergence check
         self.fp[0]        = 0.0
         self.rms          = 1.0e-40
         for i in range(0,N):
            self.fp[i]  /= self.fp[N-1]
            self.rms    += (self.fp[i] - self.fp0[i]) * (self.fp[i] - self.fp0[i])

         self.rms    = math.sqrt(self.rms)
         if it == 0:
            self.rms0   = self.rms

         # Step 9: Checking convergence
         convTest    = math.log10(self.rms) - math.log10(self.rms0)
         print("Step {:4d}: rmsLog = {:10.6f}".format(it, convTest))

         if convTest < self.tol:
            break

   # Printing results in file
   def printResults(self, outputName="resultsBlasius.dat"):
      
      N        = self.numNodes

      with open(outputName,"w") as fOut:
         fOut.write("Variables = eta, u, F, f2, fp, fpp\n")

         for i in range(0,N):
            fOut.write("{:14.10f} {:14.10f} {:14.10f} {:14.10f} {:14.10f} {:14.10f}\n".format(self.eta[i], self.fp[i] * self.VInf, self.F[i], self.f[i], self.fp[i], self.fpp[i]))

      print("Results written in file {}".format(outputName))

   def run(self):

      self.initLinear()
      self.solveBlasius()
      self.printResults()

if __name__ == '__main__':
   
      prob = Blasius(N=100,numIter=100)
      prob.run()