## example8_8
import numpy
from LUdecomp5 import *

def equations(x,h,m): # Set up finite difference eqs.
    h4 = h**4
    d = numpy.ones((m + 1),dtype = float)*6.0
    e = numpy.ones((m),dtype = float)*(-4.0)
    f = numpy.ones((m-1),dtype = float)
    b = numpy.zeros((m+1),dtype=float)
    d[0] = 1.0         
    d[1] = 7.0
    e[0] = 0.0
    f[0] = 0.0
    d[m-1] = 7.0
    d[m] = 3.0
    b[m] = 0.5*h**3
    return d,e,f,b

xStart = 0.0        # x at left end
xStop = 0.5         # x at right end
m = 20              # Number of mesh spaces
h = (xStop - xStart)/m
x = numpy.arange(xStart,xStop + h,h)
d,e,f,b = equations(x,h,m)
d,e,f = LUdecomp5(d,e,f)
y = LUsolve5(d,e,f,b)
print "\n        x              y"
for i in range(m + 1):
    print "%14.5e %14.5e" %(x[i],y[i])
raw_input("\nPress return to exit")


