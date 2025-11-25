import numpy as np

i = np.complex128(1j)

def CoordinateTrans(a,r,theta,phi):
    n = (r-i*a)*np.exp(i*phi)*np.sin(theta)
    x = np.real(n)
    y = np.imag(n)
    z = r*np.cos(theta)
    data = np.array((x,y,z))
    return data

def CartCoordInv(a,x,y,z):
    R2 = x**2+y**2+z**2
    Ra = R2-a**2
    AZ = (a**2)*(z**2)
    b3 = np.sqrt(Ra**2 + 4*AZ)
    r2 = (Ra+b3)/2
    r = np.sqrt(r2)
    theta = np.acos(z/r)
    phi = np.atan(y/x) + np.atan(a/r)
    data = np.array((r,theta,phi))
    return data

diff0 = 1/2**32

def NPtetrad(M,a,x,y,z):
    a = -a
    data0 = CartCoordInv(a,x,y,z)
    r = data0[0]
    theta = data0[1]
    phi = data0[2]
    xyz = CoordinateTrans(-a,r,theta,phi)
    xyz_r = CoordinateTrans(-a, r+diff0, theta, phi)
    xyz_theta = CoordinateTrans(-a, r, theta+diff0, phi)
    xyz_phi = CoordinateTrans(-a, r, theta, phi+diff0)
    d_dr = -(xyz - xyz_r)/diff0
    d_dtheta = -(xyz - xyz_theta)/diff0
    d_dphi = -(xyz - xyz_phi)/diff0
    nullretard = np.array([[1, 0, 0, 0], [-1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    L = np.array((0,-1,0,0))
    N = np.array((1,(1-(2*M*r)/(r**2+(a*np.cos(theta))**2))/2,0,0))
    M = np.array((i*a*np.sin(theta),i*a*np.sin(theta),1,i/(np.sin(theta))))/(np.sqrt(2)*(r-i*a*np.cos(theta)))
    Mbar = np.conj(M)
    e0 = np.zeros((4,4))*i
    e0[0,:] = L
    e0[1,:] = N
    e0[2,:] = Mbar
    e0[3,:] = M
    e0 = np.matmul(e0, nullretard)# == e0 @ nullretard
    jacobian = np.array([[1,0,0,0],[0,d_dr[0],d_dr[1],d_dr[2]],[0,d_dtheta[0],d_dtheta[1],d_dtheta[2]],[0,d_dphi[0],d_dphi[1],d_dphi[2]]])
    e0 = np.matmul(e0, jacobian)
    return e0

def NPmetric():
    g = np.zeros((4,4))
    g[0,1] = 1
    g[1,0] = 1
    g[2,3] = -1
    g[3,2] = -1
    return g

def KerrMetric(m,a,x,y,z):
    eta = NPmetric()
    e0 = NPtetrad(m,a,x,y,z)
    g = np.zeros((4,4))*i
    for u in range(0,4):
        for v in range(0,4):
            for a0 in range(0,4):
                for b0 in range(0,4):
                    g[u,v] = g[u,v] - e0[a0,u]*e0[b0,v]*eta[a0,b0]
    g = np.real(g)
    g = np.linalg.inv(g)
    return g


def Connection(m,a,x,y,z):
    # Finite difference infinitesimal
    d = 1/2**24
    dt = d
    dx = d
    dy = d
    dz = d
    # Call the Kerr Metric
    G = KerrMetric(m,a,x,y,z)
    # Compute inverse metric
    g = np.linalg.inv(G)
    # Compute metric derivatives
    Gt = KerrMetric(m,a,x,y,z)
    Gx = KerrMetric(m,a,x+dx,y,z)
    Gy = KerrMetric(m,a,x,y+dy,z)
    Gz = KerrMetric(m,a,x,y,z+dz)
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

