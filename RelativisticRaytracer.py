import random
import numpy as np
import pygame

def Lorentz(betaX,betaY,betaZ):
    beta0 = np.sqrt(betaX**2 + betaY**2 + betaZ**2)
    gamma0 = 1 / np.sqrt(1 - beta0**2)
    delta0 = gamma0**2 / (1 + gamma0)
    Lambda0 = np.zeros((4,4))
    Lambda0[0, 0] = gamma0
    Lambda0[0, 1] = gamma0 * betaX
    Lambda0[0, 2] = gamma0 * betaY
    Lambda0[0, 3] = gamma0 * betaZ
    Lambda0[1, 0] = gamma0 * betaX
    Lambda0[2, 0] = gamma0 * betaY
    Lambda0[3, 0] = gamma0 * betaZ
    Lambda0[1, 1] = 1 + (delta0 * betaX ** 2)
    Lambda0[2, 2] = 1 + (delta0 * betaY ** 2)
    Lambda0[3, 3] = 1 + (delta0 * betaZ ** 2)
    Lambda0[1, 2] = (delta0 * betaX * betaY)
    Lambda0[2, 3] = (delta0 * betaY * betaZ)
    Lambda0[3, 1] = (delta0 * betaZ * betaX)
    Lambda0[2, 1] = (delta0 * betaX * betaY)
    Lambda0[3, 2] = (delta0 * betaY * betaZ)
    Lambda0[1, 3] = (delta0 * betaZ * betaX)
    return Lambda0

def Rand(n):
    out = np.zeros((n))
    for iter in range(0,n):
        u2 = np.random.rand()
        u1 = np.random.rand()
        r = np.sqrt(-2 * np.log(u1))
        theta = 2 * np.pi * u2
        z0 = r * np.cos(theta)
        #z1 = r * np.sin(theta)
        out[iter] = z0
    return out

def V2beta(V):
    betaX = V[1] / V[0]
    betaY = V[2] / V[0]
    betaZ = V[3] / V[0]
    return betaX,betaY,betaZ

def Interval(X):
    ds2 = X[0]**2 - (X[1]**2 + X[2]**2 + X[3]**2)
    return ds2

def Vnorm(V):
    ds2 = Interval(V)
    ds = np.sqrt(ds2)
    return V/ds

def UpdatePosition(X,V):
    out = Rand(3)
    dx = out[0]
    dy = out[1]
    dz = out[2]
    dr = np.sqrt(dx**2 + dy**2 + dz**2)
    randX = np.array([0,dx,dy,dz])/dr
    betaX, betaY, betaZ = V2beta(V)
    Lambda0 = Lorentz(betaX, betaY, betaZ)
    randX0 = np.zeros((4))
    for u in range(0,4):
        for v0 in range(0,4):
            randX0[v0] += Lambda0[u,v0] * randX[u]
    ds = np.sqrt(Interval(V))
    X = X + V + randX0*ds
    return X,V

def Buffer(X,V,n):
    data = np.zeros((n-1,4))
    for tau0 in range(1,n):
        data[n-tau0-1,:] = X - V*tau0
    return data

def Step(n):
    if n>0:
        N = 1
    if n<0:
        N = -1
    if n==0:
        N = 0
    return N

def CheckLightConeTest(event0,data0):
    for tau in range(1,len(data0)):
        ds2Now = Interval(data0[tau,:] - event0)
        ds2Next = Interval(data0[tau-1, :] - event0)
        test = np.abs(Step(ds2Now) - Step(ds2Next))
        if test>0:
            return data0[tau, :]

def CheckLightCone(event0,data0):
        test = CheckLightConeTest(event0, data0)
        if test is None:
            return np.array([0,1,1,1]) * Rand(4)*len(data0)**3 - np.array([len(data0)*5,0,0,0])
        else:
            return test

def FindIndex(data,datapoint):
    n = len(data)
    ind0 = 0
    for ind in range(0,n):
        if np.sum(data[ind,:] - datapoint) == 0:
            ind0 = ind
    return ind0

def FindInterPath(event0,event1,data):
    data0 = CheckLightCone(event0, data)
    data1 = CheckLightCone(event1, data)
    n0 = FindIndex(data, data0)
    n1 = FindIndex(data, data1)
    return n0,n1

def FindSource(X,V,data):
    event0 = X
    event1 = X+V
    n0,n1 = FindInterPath(event0, event1, data)
    Source0 = data[n0:n1,:]
    return Source0

def SourceAccelerateEvent(datapoint0,datapoint1,event):
    V = datapoint1 - datapoint0
    X = (datapoint0 + datapoint1)/2
    betaX, betaY, betaZ = V2beta(V)
    Lambda0 = Lorentz(betaX, betaY, betaZ)
    LambdaInv0 = Lorentz(-betaX, -betaY, -betaZ)
    deltaX = event - X
    deltaX0 = np.zeros((4))
    for u in range(0, 4):
        for v in range(0, 4):
            deltaX0[v] += LambdaInv0[u, v] * deltaX[u]
    return deltaX0,Lambda0

def NewtonianGravity(dX):
    x = dX[1]
    y = dX[2]
    z = dX[3]
    r = np.sqrt(x**2+y**2+z**2)
    Ax = -x / r ** 3
    Ay = -y / r ** 3
    Az = -z / r ** 3
    A = np.array([0,Ax,Ay,Az])
    return A

def UpdateVelocity(X,V,data):
    datapoint = FindSource(X,V,data)
    if datapoint.any() == None:
        return X,V
    else:
        for ind0 in range(0,len(datapoint)-1):
            datapoint0 = datapoint[0+ind0]
            datapoint1 = datapoint[1+ind0]
            deltaX, Lambda0 = SourceAccelerateEvent(datapoint0,datapoint1,X)
            A_event = NewtonianGravity(deltaX)
            A0 = np.zeros((4))
            for u in range(0, 4):
                for v in range(0, 4):
                    A0[v] += Lambda0[u, v] * A_event[u]
            V = V + A0
        return X,V

def UpdateTrajectory(X,V,x,v,n0):
    Data = Buffer(X,V,int(2*n0/V[0]))
    data = Buffer(x,v,int(2*n0/v[0]))
    T = X[0]
    t = x[0]
    Tau = 0
    tau = 0
    while T<n0:
        tau = tau + 1
        Data = np.vstack((Data, X))
        X, V = UpdatePosition(X, V)
        X, V = UpdateVelocity(X, V, data)
        T = X[0]
    while t<n0:
        Tau = Tau + 1
        data = np.vstack((data, x))
        x, v = UpdatePosition(x, v)
        x, v = UpdateVelocity(x, v, Data)
        t = x[0]
    return Data,data
# or np.append(Data, X, axis=0)
X = np.array([0,0,0,0])
x = np.array([0,1000,0,0])
V = np.array([0.1,0,0,0])
v = np.array([1,0,0.1,0])
print(UpdateTrajectory(x,v,X,V,10))





















