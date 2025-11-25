import numpy as np


def Metric(M,r,theta):
    l0 = 1 - M / r
    g = np.array([[-l0,1,0,0],
                  [1,0,0,0],
                  [0,0,r**2,0],
                  [0,0,0,r**2*np.sin(theta)**2]])
    return g

def null_ray():
    p = np.array([0,-1,0,0])
    return p

def CheckRetarded(M,r,theta):
    g = Metric(M,r,theta)
    p = null_ray()
    ds2 = 0
    for mu in range(0,4):
        for nu in range(0,4):
            ds2 += g[mu,nu]*p[mu]*p[nu]
    return ds2

def InverseMetric(M,r,theta):
    g = Metric(M,r,theta)
    G = np.linalg.inv(g)
    return G

def Volume(M,r,theta):
    g = np.linalg.det(Metric(M,r,theta))
    return np.sqrt(-g)


