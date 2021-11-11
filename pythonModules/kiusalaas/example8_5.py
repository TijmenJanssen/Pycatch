## example8_5
import numpy
import run_kut5
import newtonRaphson2
import printSoln

def initCond(u):  # Initial values of [y,y',y",y"'];
                  # use 'u' if unknown
    return numpy.array([0.0, 0.0, u[0], u[1]])

def r(u):  # Boundary condition residuals-- see Eq. (8.7)
    r = numpy.zeros(len(u),dtype=float)
    X,Y = run_kut5.integrate(F,x,initCond(u),xStop,h)
    y = Y[len(Y) - 1]
    r[0] = y[2]             
    r[1] = y[3] - 1.0       
    return r
              
def F(x,y):  # First-order differential equations                   
    F = numpy.zeros((4),dtype=float)
    F[0] = y[1]
    F[1] = y[2]
    F[2] = y[3]
    if x == 0.0: F[3] = -12.0*y[1]*y[0]**2
    else:        F[3] = -4.0*(y[0]**3)/x
    return F

x = 0.0                    # Start of integration        
xStop = 1.0                # End of integration
u = numpy.array([-1.0, 1.0])     # Initial guess for u
h = 0.1                    # Initial step size
freq = 1                   # Printout frequency
u = newtonRaphson2.newtonRaphson2(r,u,1.0e-5)
X,Y = run_kut5.integrate(F,x,initCond(u),xStop,h)
printSoln.printSoln(X,Y,freq)
raw_input("\nPress return to exit")
