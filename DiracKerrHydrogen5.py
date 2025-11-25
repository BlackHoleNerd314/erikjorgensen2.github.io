import numpy as np
from scipy.linalg import expm

T0 = np.array([[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]])
X0 = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]])
Y0 = np.array([[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]])
Z0 = np.array([[0,-1,0,0],[-1,0,0,0],[0,0,0,-1],[0,0,-1,0]])

gamma = np.zeros((4,4,4))
gamma[0,:,:] = T0
gamma[1,:,:] = X0
gamma[2,:,:] = Y0
gamma[3,:,:] = Z0

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

diff0 = 1e-9

def Connection(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    # Finite difference infinitesimal
    d = diff0
    dt = d
    dx = d
    dy = d
    dz = d
    # Call the Kerr Metric
    G = KerrMetric(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    # Compute inverse metric
    g = np.linalg.inv(G)
    # Compute metric derivatives
    Gt = KerrMetric(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gx = KerrMetric(m,a,x+dx,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gy = KerrMetric(m,a,x,y+dy,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gz = KerrMetric(m,a,x,y,z+dz,Kx,Ky,Kz,Jx,Jy,Jz)
    Gt0 = KerrMetric(m, a, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gx0 = KerrMetric(m, a, x - dx, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gy0 = KerrMetric(m, a, x, y - dy, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gz0 = KerrMetric(m, a, x, y, z - dz,Kx,Ky,Kz,Jx,Jy,Jz)
    dG = np.zeros((4,4,4))
    dG[:,:,0] = (Gt-Gt0)/(-2*dt)
    dG[:,:,1] = (Gx-Gx0)/(-2*dx)
    dG[:,:,2] = (Gy-Gy0)/(-2*dy)
    dG[:,:,3] = (Gz-Gz0)/(-2*dz)
    # Corrects the finite difference ordering
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

def SpinConnect(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    omega = np.zeros((4,4,4))
    d = diff0
    dt = d
    dx = d
    dy = d
    dz = d
    GAMMA = Connection(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    e0 = Tetrad(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    omega1 = omega
    omega2 = omega
    for nu in range(0,4):
        for mu in range(0,4):
            for B in range(0,4):
                for sigma in range(0,4):
                    omega1[nu,mu,B] = omega1[nu,mu,B] + GAMMA[nu,sigma,mu]*e0[B,sigma]
    et = Tetrad(m, a, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ex = Tetrad(m, a, x + dx, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ey = Tetrad(m, a, x, y + dy, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ez = Tetrad(m, a, x, y, z + dz,Kx,Ky,Kz,Jx,Jy,Jz)
    et0 = Tetrad(m, a, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ex0 = Tetrad(m, a, x - dx, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ey0 = Tetrad(m, a, x, y - dy, z,Kx,Ky,Kz,Jx,Jy,Jz)
    ez0 = Tetrad(m, a, x, y, z - dz,Kx,Ky,Kz,Jx,Jy,Jz)
    for nu in range(0,4):
        for B in range(0,4):
            omega2[nu,0,B] = (et[B,nu]-et0[B,nu])/(2*dt)
            omega2[nu, 1, B] = (ex[B, nu] - ex0[B, nu]) / (2 * dx)
            omega2[nu, 2, B] = (ey[B, nu] - ey0[B, nu]) / (2 * dy)
            omega2[nu, 3, B] = (ez[B, nu] - ez0[B, nu]) / (2 * dz)
    omega0 = omega1+omega2
    G = KerrMetric(m, a, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    E0 = np.zeros((4,4))
    for A in range(0,4):
        for mu in range(0,4):
            for nu in range(0,4):
                E0[A,nu] = E0[A,nu] + e0[A,mu]*G[mu,nu]
    for mu in range(0,4):
        for A in range(0,4):
            for B in range(0,4):
                for nu in range(0,4):
                    omega[mu,A,B] = omega[mu,A,B] + E0[A,nu]*omega0[nu,mu,B]
    return omega

SIGMA = np.zeros((4,4,4,4))
for A in range(0,4):
    for B in range(0,4):
        SIGMA[A,B,:,:] = (np.matmul(gamma[A,:,:],gamma[B,:,:]) - np.matmul(gamma[B,:,:],gamma[A,:,:]))/4

gammaVol = T0 @ X0 @ Y0 @ Z0

def Dirac(m,a,x,y,z,d_dT,d_dX,d_dY,d_dZ,Kx,Ky,Kz,Jx,Jy,Jz):
    omega = SpinConnect(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    I = np.eye(4)
    eta = MinkowskiMetric()
    GammaConnect = np.zeros((4,4,4))
    Partial0 = np.zeros((4,4,4))
    for mu in range(0,4):
        for A in range(0,4):
            for B in range(0,4):
                for A0 in range(0,4):
                    for B0 in range(0,4):
                        GammaConnect[mu,:,:] = GammaConnect[mu,:,:] + omega[mu,A,B]*eta[A,A0]*eta[B,B0]*SIGMA[A0,B0,:,:]
    Partial0[0,:,:] = d_dT*I
    Partial0[1, :, :] = d_dX*I
    Partial0[2, :, :] = d_dY*I
    Partial0[3, :, :] = d_dZ*I
    D = Partial0 - GammaConnect
    e0 = Tetrad(m, a, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gamma = np.zeros((4,4,4))
    for A in range(0,4):
        for mu in range(0,4):
            for B in range(0,4):
                Gamma[mu,:,:] = Gamma[mu,:,:] + e0[A, mu]*eta[A,B]*gamma[B,:,:]
    D1 = np.zeros((4,4,4))
    for mu in range(0,4):
        D1[mu, :, :] = np.matmul(D[mu, :, :], Gamma[mu, :, :])
    D0 = np.zeros((4,4))
    for mu in range(0,4):
        D0[:,:] = D0[:,:] + D1[mu, :, :]
    O = np.zeros((4,4))
    O[:,:] = m*gammaVol[:,:] - D0[:,:]
    return O






# Tests
def MinkowskiMetricTest(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e0 = np.linalg.inv(Tetrad(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz))
    g = KerrMetric(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    h = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            for a in range(0,4):
                for b in range(0,4):
                    h[u,v] = h[u,v] + g[a,b]*e0[u,a]*e0[v,b]
    eta0 = MinkowskiMetric()
    print(h-eta0)
def KerrMetricTest(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    N = MinkowskiMetric()
    r = radius(a,x,y,z)
    H = Potential(m,a,r,z)
    L = Vector(a,r,x,y,z)
    # Construct Kerr Spacetime
    G = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            G[u,v] = N[u,v] + H*L[v]*L[u]
    G0 = KerrMetric(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    print(G-G0)
def GammaFrameTest():
    eta = MinkowskiMetric()
    I = np.eye(4)
    for mu in range(0,4):
        for nu in range(0,4):
            Amat = gamma[mu,:,:]
            Bmat = gamma[nu,:,:]
            Cmat = np.matmul(Amat,Bmat)
            Cmat0 = np.matmul(Bmat,Amat)
            mat = (Cmat + Cmat0)/2
            Imat0 = I * eta[mu,nu]
            print(mat-Imat0)
#GammaFrameTest()
m = 1.3
a = 0.13
x = 3
y = 4
z = 12
Kx = 1
Ky = 2
Kz = 3
Jx = 4
Jy = 5
Jz = 6
KerrMetricTest(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
#print(test) # Should be np.zeros((4,4))
MinkowskiMetricTest(m,a,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
#print(test0) # Should be np.zeros((4,4))

#e0 = Tetrad(m, a, x, y, z, Kx, Ky, Kz, Jx, Jy, Jz)
#print("Determinant of e0:", np.linalg.det(e0))  # Should be 1 (for the tetrad)
#g = KerrMetric(m, a, x, y, z, Kx, Ky, Kz, Jx, Jy, Jz)
#print("Determinant of g:", np.linalg.det(g))    # Should be -1 (det(η) in 4D)

print(SpinConnect(13,1,3,4,12, Kx, Ky, Kz, Jx, Jy, Jz))
x = 40
y = 30
z = 120
d_dt = 0
d_dx = 0
d_dy = 0
d_dz = 0
M = 1
m = M/2
a = 1/(2*M)
Operator = Dirac(m,a,x,y,z,d_dt,d_dx,d_dy,d_dz,Kx,Ky,Kz,Jx,Jy,Jz)
print(Operator)