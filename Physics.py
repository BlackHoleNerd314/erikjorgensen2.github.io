import numpy as np
from CollapseFeynman import Factorial0

def P(m,l,theta):
    term0 = (-1)**m
    term1 = (np.sin(theta))**m
    term2 = term0 * term1
    term12 = 0
    for k in range(0,np.floor(l/2)):
        term3 = (-1)**k
        term4 = Factorial0(2*l-2*k)
        term5 = (2**l)*Factorial0(k)
        term6 = Factorial0(l-k)
        term9 = (l - 2 * k - m)
        term7 = Factorial0(term9)
        term8 = term4/(term5*term6*term7)
        term10 = (np.cos(theta))**term9
        term11 = term3 * term8 * term10
        term12 = term12 + term11
    term13 = term2 * term12
    return term13

def L(alpha,n,x):
    term9 = 0
    for k in range(0,n):
        term1 = (-1)**k
        term2 = Factorial0(n+alpha)
        term3 = Factorial0(n-k)
        term4 = Factorial0(alpha+k)
        term5 = Factorial0(k)
        term6 = term2/(term3*term4*term5)
        term7 = x**k
        term8 = term1 * term6 * term7
        term9 = term9 + term8
    return term9