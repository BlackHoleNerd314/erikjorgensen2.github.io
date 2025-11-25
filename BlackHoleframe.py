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

def Tetrad(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    X0 = np.array([t,x,y,z])
    delta0 = LorentzTransform(Kx,Ky,Kz,Jx,Jy,Jz)
    X = np.linalg.inv(delta0) @ X0
    eta = MinkowskiMetric()
    x0 = X[1]
    y0 = X[2]
    z0 = X[3]
    r = radius(a,x0,y0,z0)
    H = Potential(m,a,r,z0)
    L = Vector(a,r,x0,y0,z0)
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

def KerrMetric(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e0 = Tetrad(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    eta = MinkowskiMetric()
    g = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            for a in range(0,4):
                for b in range(0,4):
                    g[u,v] = g[u,v] + eta[a,b]*e0[u,a]*e0[v,b]
    return g

diff0 = 1e-9

def Connection(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    # Finite difference infinitesimal
    d = diff0
    dt = d
    dx = d
    dy = d
    dz = d
    # Call the Kerr Metric
    G = KerrMetric(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    # Compute inverse metric
    g = np.linalg.inv(G)
    # Compute metric derivatives
    Gt = KerrMetric(m,a,t+dt,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gx = KerrMetric(m,a,t,x+dx,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gy = KerrMetric(m,a,t,x,y+dy,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gz = KerrMetric(m,a,t,x,y,z+dz,Kx,Ky,Kz,Jx,Jy,Jz)
    Gt0 = KerrMetric(m, a,t - dt, x, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gx0 = KerrMetric(m, a,t, x - dx, y, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gy0 = KerrMetric(m, a,t, x, y - dy, z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gz0 = KerrMetric(m, a,t, x, y, z - dz,Kx,Ky,Kz,Jx,Jy,Jz)
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

# WIP

def spin_connection(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    # e_mu_a: tetrad with first index mu (row) and second index a (col) as in your code
    e = Tetrad(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    inv_e = np.linalg.inv(e)   # inverse tetrad
    Gam = Connection(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)  # Gamma^o_{uv} as provided
    # compute partial_mu e^lambda{}_b using finite difference (reuse diff0)
    d = diff0
    dE = np.zeros((4,4,4))  # dE[mu,lambda,b]
    for mu in range(0,4):
        for u in range(0,4):
            for a0 in range(0,4):
                Xp = [t,x,y,z].copy()
                Xp[mu] += d
                Ep = Tetrad(m,a,Xp[0],Xp[1],Xp[2],Xp[3],Kx,Ky,Kz,Jx,Jy,Jz)
                Xm = [t,x,y,z].copy()
                Xm[mu] -= d
                Em = Tetrad(m,a,Xm[0],Xm[1],Xm[2],Xm[3],Kx,Ky,Kz,Jx,Jy,Jz)
                dE[mu,u,a0] = (Ep[u,a0] - Em[u,a0])/(2*d)
    # now compute omega_mu^{ab}
    omega = np.zeros((4,4,4))
    eta = MinkowskiMetric()
    G = KerrMetric(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    term1 = dE #np.zeros((4,4,4))
    for mu in range(0,4):
        for nu in range(0,4):
            for B in range(0,4):
                for sigma in range(0,4):
                    term1[mu,nu,B] += + Gam[nu,sigma,mu] * e[sigma,B]
    e1 = np.zeros((4,4))
    for A in range(0,4):
        for u in range(0,4):
            for v in range(0,4):
                e1[v,A] += e[u,A] * G[u,v]
    for A in range(0,4):
        for B in range(0,4):
            for mu in range(0,4):
                for nu in range(0,4):
                    omega[mu,A,B] += e1[nu,A] * term1[mu,nu,B]
    # impose antisymmetry in a,b if needed: (omega_{mu}^{ab} = -omega_{mu}^{ba})
    return omega
### wrong
def curved_gammas(e):
    inv_e = np.linalg.inv(e)   # e^mu_a
    gammas_frame = gamma  # your gamma[0..3,:,:] arrays
    G = np.zeros((4,4,4), dtype=float)
    for mu in range(4):
        mat = np.zeros((4,4))
        for a_idx in range(4):
            mat += inv_e[mu,a_idx] * gammas_frame[a_idx,:,:]
        G[mu,:,:] = mat
    return G

def spinor_connection_matrices(omega):
    # returns S_mu (4 matrices 4x4)
    S = [np.zeros((4,4)) for _ in range(4)]
    for mu in range(4):
        Smu = np.zeros((4,4))
        for a_idx in range(4):
            for b_idx in range(4):
                Smu += 0.25 * omega[mu,a_idx,b_idx] * (gamma[a_idx,:,:] @ gamma[b_idx,:,:])
        S[mu] = Smu
    return S

def Dirac_operator_action(psi_vec, m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e = Tetrad(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    Gammas = curved_gammas(e)
    omega = spin_connection(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    S = spinor_connection_matrices(omega)
    # compute partial_mu psi (finite differences or spectral)
    # Here show symbolic local action: (i * Gamma^mu (partial_mu + S_mu) - m) psi
    # Implement finite diff for time evolution separately (see next)
    #return #operator_matrix  # you can build local 4x4 matrix per derivative + multiply by derivatives

##

print(spin_connection(0.5,0.1,0,10,1,0.1,0,0,0,0,0,0))

##

import numpy as np

def asym_norm(omega):
    A = omega + np.swapaxes(omega,1,2)
    return np.max(np.abs(A)), np.argmax(np.abs(A)), A

for d in [1e-8, 5e-9, 2.5e-9, 1.25e-9]:
    diff0 = d
    # make sure your spin_connection uses the global diff0 or pass it in
    omega = spin_connection(m=1.0, a=0.5,
                            t=0.0, x=2.0, y=0.1, z=0.2,
                            Kx=0,Ky=0,Kz=0,Jx=0,Jy=0,Jz=0)
    mx, idx, A = asym_norm(omega)
    print("d=", d, "max asym=", mx, "idx(flat)=", idx)

# also print per-mu maxima
for mu in range(4):
    rowmax = np.max(np.abs(omega[mu] + omega[mu].T))
    print("mu",mu,"asym max",rowmax)






# Tests
def MinkowskiMetricTest(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e0 = np.linalg.inv(Tetrad(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz))
    g = KerrMetric(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
    h = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            for a in range(0,4):
                for b in range(0,4):
                    h[u,v] = h[u,v] + g[a,b]*e0[u,a]*e0[v,b]
    eta0 = MinkowskiMetric()
    print(h-eta0)
def KerrMetricTest(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    X0 = np.array([t, x, y, z])
    delta0 = LorentzTransform(Kx, Ky, Kz, Jx, Jy, Jz)
    X = np.linalg.inv(delta0) @ X0
    x0 = X[1]
    y0 = X[2]
    z0 = X[3]
    N = MinkowskiMetric()
    r = radius(a,x0,y0,z0)
    H = Potential(m,a,r,z0)
    L = Vector(a,r,x0,y0,z0)
    # Construct Kerr Spacetime
    G = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            G[u,v] = N[u,v] + H*L[v]*L[u]
    G0 = KerrMetric(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz)
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
def VolumeTest(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    e0 = Tetrad(m, a, t, x, y, z, Kx, Ky, Kz, Jx, Jy, Jz)
    G = KerrMetric(m, a, t, x, y, z, Kx, Ky, Kz, Jx, Jy, Jz)
    print(np.linalg.det(e0)-1)
    print(np.linalg.det(G)+1)
    print(np.linalg.det(T0)-1)
    print(np.linalg.det(X0)-1)
    print(np.linalg.det(Y0)-1)
    print(np.linalg.det(Z0)-1)
def LorentzTetradTest(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    G0 = (Tetrad(m,a,t,x,y,z,0,0,0,0,0,0))
    c0 = (LorentzTransform(Kx,Ky,Kz,Jx,Jy,Jz))
    X = c0 @ np.array([t,x,y,z])#@ np.linalg.inv(c0)
    cG0 = Tetrad(m,a,X[0],X[1],X[2],X[3],Kx,Ky,Kz,Jx,Jy,Jz)
    cG1 = G0 @ c0
    print(cG0-cG1)


