## example8_4
import numpy
from  bulStoer import *
import newtonRaphson2
import printSoln

def initCond(u):  # Initial values of [y,y',y",y"'];
                  # use 'u' if unknown
    return numpy.array([0.0, u[0], 0.0, u[1]])

def r(u):  # Boundary condition residuals--see Eq. (8.7)
    r = numpy.zeros(len(u),dtype=float)
    X,Y = bulStoer(F,xStart,initCond(u),xStop,H)
    y = Y[len(Y) - 1]
    r[0] = y[0]   
    r[1] = y[2]
    return r
              
def F(x,y):  # First-order differential equations                   
    F = numpy.zeros((4),dtype=float)
    F[0] = y[1]
    F[1] = y[2]
    F[2] = y[3]
    F[3] = x
    return F

xStart = 0.0               # Start of integration        
xStop = 1.0                # End of integration
u = numpy.array([0.0, 1.0])      # Initial guess for {u}
H = 0.5                    # Printout incremant
freq = 1                   # Printout frequency
u = newtonRaphson2.newtonRaphson2(r,u,1.0e-4)
X,Y = bulStoer(F,xStart,initCond(u),xStop,H)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
