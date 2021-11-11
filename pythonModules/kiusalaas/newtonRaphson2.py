## module newtonRaphson2
''' soln = newtonRaphson2(f,x,tol=1.0e-9).
    Solves the simultaneous equations f(x) = 0 by
    the Newton-Raphson method using {x} as the initial
    guess. Note that {f} and {x} are vectors.
'''
import numpy
import math
import gaussPivot

def newtonRaphson2(f,x,tol=1.0e-9):
    
    def jacobian(f,x):
        h = 1.0e-4
        n = len(x)
        jac = numpy.zeros((n,n),dtype=float)
        f0 = f(x)
        for i in range(n):
            temp = x[i]
            x[i] = temp + h
            f1 = f(x)
            x[i] = temp
            jac[:,i] = (f1 - f0)/h
        return jac,f0
    
    for i in range(30):
        jac,f0 = jacobian(f,x)
        if math.sqrt(numpy.dot(f0,f0)/len(x)) < tol: return x
        dx = gaussPivot.gaussPivot(jac,-f0)
        x = x + dx
        if math.sqrt(numpy.dot(dx,dx)) < tol*max(max(abs(x)),1.0): return x
    print 'Too many iterations'
