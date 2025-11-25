import numpy as np
import math
from scipy.linalg import expm

def LorentzTransform(Kx,Ky,Kz,Jx,Jy,Jz):
    # Define general Lorentz Transformation matrix
    M = np.array([[0,Kx,Ky,Kz],
                  [Kx,0,-Jz,Jy],
                  [Ky,Jz,0,-Jx],
                  [Kz,-Jy,Jx,0]])
    # Compute Matrix Exponential of Generating Matrix
    Lorentz = expm(M)
    return Lorentz

def MinkowskiMetric():
    # The spacetime invariant interval
    N = np.zeros((4,4))
    N[0,0] = -1
    N[1,1] = 1
    N[2,2] = 1
    N[3,3] = 1
    return N

def radius(a,x,y,z):
    R2 = x*x + y*y + z*z
    Ra = R2 - a*a
    AZ = a*a*z*z
    b3 = np.sqrt(Ra*Ra + 4*AZ)
    r2 = (Ra + b3)/2
    r = np.sqrt(r2)
    return r

def Potential(m,a,r,z):
    AZ = a * a * z * z
    H1 = 2 * m * r * r * r
    H2 = r * r * r * r + AZ
    H = H1 / H2
    return H

def Vector(a,r,x,y,z):
    L = np.zeros((4, 1))
    L[0] = 1
    L[1] = (r * x + a * y) / (r * r + a * a)
    L[2] = (r * y - a * x) / (r * r + a * a)
    L[3] = z / r
    return L

def Tetrad(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    delta0 = LorentzTransform(Kx,Ky,Kz,Jx,Jy,Jz)
    eta = MinkowskiMetric()
    r = radius(a,x,y,z)
    H = Potential(m,a,r,z)
    L = Vector(a,r,x,y,z)
    e0 = delta0
    l = np.zeros((4))
    for v in range(0,4):
        for a in range(0, 4):
            l[v] = l[v] + eta[a,v]*L[a]
    L0 = np.zeros((4))
    for n1 in range(0,4):
        for n2 in range(0,4):
            L0[n1] = L0[n1] + l[n2] * delta0[n2,n1]
    for u in range(0,4):
        for v in range(0,4):
            e0[u,v] = e0[u,v] + L0[v]*L[u]*H/2

    return e0

def KerrMetric(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e0 = Tetrad(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    eta = MinkowskiMetric()
    g = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            for a in range(0,4):
                for b in range(0,4):
                    g[u,v] = g[u,v] + eta[a,b]*e0[u,a]*e0[v,b]
    return g

def IntervalFlat(X,eta):
    eta = MinkowskiMetric()
    s2 = 0
    for u in range(0,4):
        for v in range(0,4):
            s2 = s2 + eta[u,v]*X[u]*X[v]
    s = np.sqrt(-s2)
    return s

def FlatNormalized(X):
    eta = MinkowskiMetric()
    s = IntervalFlat(X,eta)
    X = X/s
    return X

def PlaneWave(X,P):
    eta = MinkowskiMetric()
    phase = 0
    pi = math.pi
    for u in range(0,4):
        for v in range(0,4):
            phase = phase + eta[u,v]*X[v]*P[u]
    cycles = (2j*pi*phase)
    return cycles

def NewtonGravity(m,X):
    r = np.sqrt(X[1] ** 2 + X[2] ** 2 + X[3] ** 2)
    A = m * X / r ** 3
    A[0] = 0
    return A

def Connection(m,a,x,y,z):
    # Finite difference infinitesimal
    d = 1/2**24
    dt = d
    dx = d
    dy = d
    dz = d
    # Call the Kerr Metric
    G = KerrMetric(m,a,x,y,z,0,0,0,0,0,0)
    # Compute inverse metric
    g = np.linalg.inv(G)
    # Compute metric derivatives
    Gt = KerrMetric(m,a,x,y,z,0,0,0,0,0,0)
    Gx = KerrMetric(m,a,x+dx,y,z,0,0,0,0,0,0)
    Gy = KerrMetric(m,a,x,y+dy,z,0,0,0,0,0,0)
    Gz = KerrMetric(m,a,x,y,z+dz,0,0,0,0,0,0)
    dG = np.zeros((4,4,4))
    dG[:,:,0] = (Gt-G)/dt
    dG[:,:,1] = (Gx-G)/dx
    dG[:,:,2] = (Gy-G)/dy
    dG[:,:,3] = (Gz-G)/dz
    # Corrects the finite difference ordering
    dG = - dG
    B1 = np.zeros((4,4,4))
    B2 = np.zeros((4, 4, 4))
    B3 = np.zeros((4, 4, 4))
    # Components of the Connection using metric derivatives
    for u in range(0,4):
        for v in range(0,4):
            for p in range(0,4):
                B1[u,v,p]=dG[v,p,u]
                B2[u,v,p]=dG[p,u,v]
                B3[u,v,p]=-dG[u,v,p]
    B = B1+B2+B3
    # Compute the Covariant Derivative Connection Coefficients
    L12 = np.zeros((4,4,4))
    for o in range(0,4):
        for p in range(0,4):
            for u in range(0,4):
                for v in range(0,4):
                    L12[o,u,v]=L12[o,u,v]+(g[o,p]*B[u,v,p]/2)
    return L12

def Geodesic(m,a,X,V,J,n):
    data = np.zeros((n,12))
    # Initialize Trajectory and Spin Axis
    for u in range(0,4):
        data[0,u] = X[u]
        data[0,u+4] = V[u]
        data[0,u+8] = J[u]
    # Iterate the trajectory in proper time
    for s in range(1,n):
        # gravitating mass is stationary, so we call the coordinates directly
        x = X[1]
        y = X[2]
        z = X[3]
        # Call the Metric (Geometry) and Connection (Forces)
        G = KerrMetric(m,a,x,y,z,0,0,0,0,0,0)
        L = Connection(m,a,x,y,z)
        # Normalize the mass and 4-momentum to compute 4-velocity
        m2 = 0
        for u in range(0,4):
            for v in range(0,4):
                m2 += V[u]*G[u,v]*V[v]
        V = V/np.sqrt(-m2)
        # Compute parallel transport of 4-velocity and 4-spin vectors
        A = np.zeros((4))
        dJ = np.zeros((4))
        for u in range(0,4):
            for v in range(0,4):
                for o in range(0,4):
                    A[o] += L[o,u,v]*V[u]*V[v]
                    dJ[o] += L[o,u,v]*V[u]*J[v]
        # Update all 4-vectors
        V = V+A
        J = J+dJ
        X = X+V
        # Update input parameters
        for p in range(0,4):
            data[s,p] = X[p]
            data[s,p+4] = V[p]
            data[s,p+8] = J[p]
    return data

def o(X,V):
    V[0] = 1
    X = X + V
    return X,V

def C(X,V):
    V = FlatNormalized(V)
    X = X + V
    return X,V

def G(m,X,V):
    V = V - NewtonGravity(m, X)
    V[0] = 1
    X = X + V
    return m,X,V

def CG(m,X,V):
    data = Geodesic(m,0,X,V,np.array([0,0,0,0]),2)
    X = data[1,0:4]
    V = data[1,4:8]
    return m,X,V

def H(m,X,V):
    V[0] = 1
    P = m*V
    normP = np.sqrt(P[1] ** 2 + P[2] ** 2 + P[3] ** 2)
    P[0] = normP**2 / (2*m)
    phase = PlaneWave(X,P)
    Psi = np.exp(phase)
    return Psi

def CH(m,X,V):
    V = FlatNormalized(V)
    P = m*V
    phase = PlaneWave(X, P)
    Psi = np.exp(phase)
    return Psi

#def GH():
print(CG(1,np.array([0,100,0,0]),np.array([1,0,0.1,0])))
