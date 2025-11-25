import math
import numpy as np
from scipy.linalg import expm

def MinkowskiMetric():
    # The spacetime invariant interval
    N = np.zeros((4,4))
    N[0,0] = -1
    N[1,1] = 1
    N[2,2] = 1
    N[3,3] = 1
    return N

def KerrMetric(m,a,x,y,z):
    # Standard Minkowski Spacetime
    N = MinkowskiMetric()
    # Define Radius
    R2 = x*x + y*y + z*z
    Ra = R2 - a*a
    AZ = a*a*z*z
    b3 = np.sqrt(Ra*Ra + 4*AZ)
    r2 = (Ra + b3)/2
    r = np.sqrt(r2)
    # Define Scalar Perturbation
    H1 = 2*m*r*r*r
    H2 = r*r*r*r + AZ
    H = H1/H2
    # Define Vector Perturbation
    Lt = 1
    Lx = (r*x + a*y)/(r*r+a*a)
    Ly = (r*y - a*x)/(r*r+a*a)
    Lz = z/r
    L = np.zeros((4,1))
    L[1] = Lx
    L[2] = Ly
    L[3] = Lz
    L[0] = Lt
    # Construct Kerr Spacetime
    G = np.zeros((4,4))
    for u in range(0,4):
        for v in range(0,4):
            G[u,v] = N[u,v] + H*L[v]*L[u]
    return G

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

def Geodesic(m,a,X,V,J,n):
    data = np.zeros((n,8))
    # Initialize Trajectory and Spin Axis
    for u in range(0,4):
        data[0,u] = X[u]
        data[0,u+4] = J[u]
    # Iterate the trajectory in proper time
    for s in range(1,n):
        # gravitating mass is stationary, so we call the coordinates directly
        x = X[1]
        y = X[2]
        z = X[3]
        # Call the Metric (Geometry) and Connection (Forces)
        G = KerrMetric(m,a,x,y,z)
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
            data[s,p+4] = J[p]
    return data

def LorentzTransform(Kx,Ky,Kz,Jx,Jy,Jz):
    # Define general Lorentz Transformation matrix
    M = np.array([[0,Kx,Ky,Kz],
                  [Kx,0,-Jz,Jy],
                  [Ky,Jz,0,-Jx],
                  [Kz,-Jy,Jx,0]])
    # Compute Matrix Exponential of Generating Matrix
    Lorentz = expm(M)
    return Lorentz

def Boost(V):
    # Normalize 4-Velocity
    G = MinkowskiMetric()
    M2 = 0
    for u in range(0,4):
        for v in range(0,4):
            M2 += V[u] * G[u,v] * V[v]
    V = V / np.sqrt(-M2)
    # Convert from 3-velocity to Rapidity (Useful for Boosts)
    K0_norm = np.sqrt(V[1] ** 2 + V[2] ** 2 + V[3] ** 2)
    boostT = np.zeros((4,4))
    # Error checking for identity
    if K0_norm == 0:
        boostT = np.eye(4,4)
    # Lorentz Boost matrix for the 4-velocity in rest frame
    if K0_norm != 0:
        Kx = (np.acosh(V[0])) * V[1] / K0_norm
        Ky = (np.acosh(V[0])) * V[2] / K0_norm
        Kz = (np.acosh(V[0])) * V[3] / K0_norm
        boostT = LorentzTransform(-Kx, -Ky, -Kz, 0, 0, 0)
    return boostT

def Rotate(A):
    # Error checking for identity
    if A[1]**2 + A[2]**2 == 0:
        rotZ = np.eye(4,4)
    else:
    # Compute necessary normalized cross products
        A_norm = np.sqrt(A[1]**2+A[2]**2+A[3]**2)
        A = A / A_norm
        Jx = A[2]*1 - A[3]*0
        Jy = A[3]*0 - A[1]*1
        Jz = A[1]*0 - A[2]*0
        J_norm = np.sqrt(Jx ** 2 + Jy ** 2 + Jz ** 2)
    # Error correct for negative angles in euler Z axis
        if A[3] < 0:
            Jx0 = (math.pi - np.asin(J_norm)) * Jx / J_norm
            Jy0 = (math.pi - np.asin(J_norm)) * Jy / J_norm
            Jz0 = (math.pi - np.asin(J_norm)) * Jz / J_norm
    # Convert to Euler angles (Useful for Rotations)
        else:
            Jx0 = np.asin(J_norm) * Jx / J_norm
            Jy0 = np.asin(J_norm) * Jy / J_norm
            Jz0 = np.asin(J_norm) * Jz / J_norm
        rotZ = LorentzTransform(0,0,0, -Jx0, -Jy0, -Jz0)
    return rotZ

def LightCone(m,a,X,V,J,n):
    # Compute Geodesic (Orbital Trajectory)
    orbit = Geodesic(m,a,X,V,J,n)
    # Initialize Reference Proper Time
    s = 0
    t = orbit[s, 0]
    x = orbit[s, 1]
    y = orbit[s, 2]
    z = orbit[s, 3]
    # Initialize Kerr Radial Coordinate
    R2 = x * x + y * y + z * z
    Ra = R2 - a * a
    AZ = a * a * z * z
    b3 = np.sqrt(Ra * Ra + 4 * AZ)
    r2 = (Ra + b3) / 2
    r = np.sqrt(r2)
    # Kerr Radial Coordinate defines Ingoing Geodesic Time
    while t<r:
        s += 1
        t = orbit[s,0]
        x = orbit[s,1]
        y = orbit[s,2]
        z = orbit[s,3]
        # Define Kerr Radial Coordinate per Interation
        R2 = x * x + y * y + z * z
        Ra = R2 - a * a
        AZ = a * a * z * z
        b3 = np.sqrt(Ra * Ra + 4 * AZ)
        r2 = (Ra + b3) / 2
        r = np.sqrt(r2)
    return s

def RetardedTime(m,a,X,V,J,n):
    # Call the Geodesic Calculator
    orbit = Geodesic(m, a, X, V, J, n)
    s = LightCone(m, a, X, V, J, n)
    # Use the Geodesics as a Parameterized Trajectory
    event0 = orbit[s]
    t = event0[0]
    x = event0[1]
    y = event0[2]
    z = event0[3]
    Jt = event0[4]
    Jx = event0[5]
    Jy = event0[6]
    Jz = event0[7]
    # Compute 4-velocity from 4-Position
    event1 = orbit[s-1]
    Pt = event0[0]-event1[0]
    Px = event0[1]-event1[1]
    Py = event0[2]-event1[2]
    Pz = event0[3]-event1[3]
    # Return Event Data on Ingoing Light Cone
    event0 = np.array((t,x,y,z,Pt,Px,Py,Pz,Jt,Jx,Jy,Jz))
    return event0


def FrameHop(event0,ReferenceEvent):
    t0 = event0[0]
    x0 = event0[1]
    y0 = event0[2]
    z0 = event0[3]
    Pt0 = event0[4]
    Px0 = event0[5]
    Py0 = event0[6]
    Pz0 = event0[7]
    Jt0 = event0[8]
    Jx0 = event0[9]
    Jy0 = event0[10]
    Jz0 = event0[11]
    t = ReferenceEvent[0]
    x = ReferenceEvent[1]
    y = ReferenceEvent[2]
    z = ReferenceEvent[3]
    Pt = ReferenceEvent[4]
    Px = ReferenceEvent[5]
    Py = ReferenceEvent[6]
    Pz = ReferenceEvent[7]
    Jt = ReferenceEvent[8]
    Jx = ReferenceEvent[9]
    Jy = ReferenceEvent[10]
    Jz = ReferenceEvent[11]
    t = t - t0
    x = x - x0
    y = y - y0
    z = z - z0
    Mboost = Boost(np.array((Pt0,Px0,Py0,Pz0)))
    X = np.array((t,x,y,z))
    X1 = np.zeros(4)
    P = np.array((Pt,Px,Py,Pz))
    P1 = np.zeros(4)
    J = np.array((Jt,Jx,Jy,Jz))
    J1 = np.zeros(4)
    J0 = np.array((Jt0, Jx0, Jy0, Jz0))
    J2 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += Mboost[u,v]*X[v]
            P1[u] += Mboost[u,v]*P[v]
            J1[u] += Mboost[u,v]*J[v]
            J2[u] += Mboost[u,v]*J0[v]
    X = X1
    P = P1
    J = J1
    Mrotate = np.linalg.inv(Rotate(J2))
    X1 = np.zeros(4)
    P1 = np.zeros(4)
    J1 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += Mrotate[u,v]*X[v]
            P1[u] += Mrotate[u,v]*P[v]
            J1[u] += Mrotate[u,v]*J[v]
    #M1 = np.linalg.inv(Mrotate)
    #M2 = np.linalg.inv(Mboost)
    #M3 = Mrotate
    #J = J1
    #J1 = np.zeros(4)
    #for u in range(0, 4):
    #    for v in range(0, 4):
    #        J1[u] += M1[u, v] * J[v]
    #J = J1
    #J1 = np.zeros(4)
    #for u in range(0, 4):
    #    for v in range(0, 4):
    #        J1[u] += M2[u, v] * J[v]
    #J = J1
    #J1 = np.zeros(4)
    #for u in range(0, 4):
    #    for v in range(0, 4):
    #        J1[u] += M3[u, v] * J[v]
    event1 = np.array((X1,P1,J1))
    return event1

def FrameHopInverse(event0,event1):
    t0 = event0[0]
    x0 = event0[1]
    y0 = event0[2]
    z0 = event0[3]
    Pt0 = event0[4]
    Px0 = event0[5]
    Py0 = event0[6]
    Pz0 = event0[7]
    Jt0 = event0[8]
    Jx0 = event0[9]
    Jy0 = event0[10]
    Jz0 = event0[11]

    Mboost0 = Boost(np.array((Pt0, Px0, Py0, Pz0)))
    J0 = np.array((Jt0, Jx0, Jy0, Jz0))
    J2 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            J2[u] += Mboost0[u,v]*J0[v]






    X = np.zeros(4)
    X[0] = event1[0]
    X[1] = event1[1]
    X[2] = event1[2]
    X[3] = event1[3]
    P = np.zeros(4)
    P[0] = event1[4]
    P[1] = event1[5]
    P[2] = event1[6]
    P[3] = event1[7]
    J = np.zeros(4)
    J[0] = event1[8]
    J[1] = event1[9]
    J[2] = event1[10]
    J[3] = event1[11]

    Mrotate = Rotate(J2)
    X1 = np.zeros(4)
    P1 = np.zeros(4)
    J1 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += Mrotate[u,v]*X[v]
            P1[u] += Mrotate[u,v]*P[v]
            J1[u] += Mrotate[u,v]*J[v]

    Mboost = np.linalg.inv(Boost(np.array((Pt0,Px0,Py0,Pz0))))
    X = X1
    X1 = np.zeros(4)
    P = P1
    P1 = np.zeros(4)
    J = J1
    J1 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += Mboost[u,v]*X[v]
            P1[u] += Mboost[u,v]*P[v]
            J1[u] += Mboost[u,v]*J[v]
    X = X1
    P = P1
    J = J1

    t = X[0]
    x = X[1]
    y = X[2]
    z = X[3]
    Pt = P[0]
    Px = P[1]
    Py = P[2]
    Pz = P[3]
    Jt = J[0]
    Jx = J[1]
    Jy = J[2]
    Jz = J[3]


    t = t + t0
    x = x + x0
    y = y + y0
    z = z + z0

    ReferenceEvent = np.array((t,x,y,z,Pt,Px,Py,Pz,Jt,Jx,Jy,Jz))

    return ReferenceEvent

def OrbitPair(X,V,J,x,v,j):
    m2 = 0
    a2 = 0
    M2 = 0
    A2 = 0
    eta0 = MinkowskiMetric()
    for u0 in range(0, 4):
        for v0 in range(0, 4):
            m2 += -eta0[u0, v0] * v[u0] * v[v0]
            a2 += eta0[u0, v0] * j[u0] * j[v0]
            M2 += -eta0[u0, v0] * V[u0] * V[v0]
            A2 += eta0[u0, v0] * J[u0] * J[v0]
    M = np.sqrt(M2)
    A = np.sqrt(A2)
    m = np.sqrt(m2)
    a = np.sqrt(a2)
    particle1 = np.array((x[0],x[1],x[2],x[3],v[0],v[1],v[2],v[3],j[0],j[1],j[2],j[3]))
    particle2 = np.array((X[0],X[1],X[2],X[3],V[0],V[1],V[2],V[3],J[0],J[1],J[2],J[3]))

    particle2frame1 = FrameHop(particle1,particle2)
    Xnew = np.zeros(4)
    Xnew[0] = particle2frame1[0,0]
    Xnew[1] = particle2frame1[0,1]
    Xnew[2] = particle2frame1[0,2]
    Xnew[3] = particle2frame1[0,3]
    Vnew = np.zeros(4)
    Vnew[0] = particle2frame1[1,0]
    Vnew[1] = particle2frame1[1,1]
    Vnew[2] = particle2frame1[1,2]
    Vnew[3] = particle2frame1[1,3]
    Jnew = np.zeros(4)
    Jnew[0] = particle2frame1[2,0]
    Jnew[1] = particle2frame1[2,1]
    Jnew[2] = particle2frame1[2,2]
    Jnew[3] = particle2frame1[2,3]

    particle1frame2 = FrameHop(particle2, particle1)
    xnew = np.zeros(4)
    xnew[0] = particle1frame2[0,0]
    xnew[1] = particle1frame2[0,1]
    xnew[2] = particle1frame2[0,2]
    xnew[3] = particle1frame2[0,3]
    vnew = np.zeros(4)
    vnew[0] = particle1frame2[1,0]
    vnew[1] = particle1frame2[1,1]
    vnew[2] = particle1frame2[1,2]
    vnew[3] = particle1frame2[1,3]
    jnew = np.zeros(4)
    jnew[0] = particle1frame2[2,0]
    jnew[1] = particle1frame2[2,1]
    jnew[2] = particle1frame2[2,2]
    jnew[3] = particle1frame2[2,3]

    gamma_factor = 1/np.sqrt(1-((Vnew[1] ** 2 + Vnew[2] ** 2 + Vnew[3] ** 2)/(Vnew[0]**2)))
    Gamma_Factor = 1/np.sqrt(1-((vnew[1] ** 2 + vnew[2] ** 2 + vnew[3] ** 2)/(vnew[0]**2)))
    n = (gamma_factor+Gamma_Factor) * np.sqrt(Xnew[1] ** 2 + Xnew[2] ** 2 + Xnew[3] ** 2)
    N = (gamma_factor+Gamma_Factor) * np.sqrt(xnew[1] ** 2 + xnew[2] ** 2 + xnew[3] ** 2)

    particle2frame1New = RetardedTime(m, a, Xnew, Vnew, Jnew, math.floor(n))
    particle1frame2New = RetardedTime(M, A, xnew, vnew, jnew, math.floor(N))

    #particle1New = FrameHop(particle1frame2New, particle2)
    #particle2New = FrameHop(particle2frame1New, particle1)
    particle1New = FrameHopInverse(particle2,particle1frame2New)
    particle2New = FrameHopInverse(particle1,particle2frame1New)

    particle1 = particle1New
    particle2 = particle2New

    data0 = np.array((particle2,particle1))
    return data0



    #particle1frame1 = [0,0,0,0,1,0,0,0,0,0,0,1]
    #particle2frame1 =
    #particle1frame2 = FrameHop(particle2frame1,particle1frame1)


def LinearPair(X,V,J,x,v,j):
    m = 0
    a = 0
    M = 0
    A = 0
    eta0 = MinkowskiMetric()
    particle1 = np.array((x[0],x[1],x[2],x[3],v[0],v[1],v[2],v[3],j[0],j[1],j[2],j[3]))
    particle2 = np.array((X[0],X[1],X[2],X[3],V[0],V[1],V[2],V[3],J[0],J[1],J[2],J[3]))

    particle2frame1 = FrameHop(particle1,particle2)
    Xnew = np.zeros(4)
    Xnew[0] = particle2frame1[0,0]
    Xnew[1] = particle2frame1[0,1]
    Xnew[2] = particle2frame1[0,2]
    Xnew[3] = particle2frame1[0,3]
    Vnew = np.zeros(4)
    Vnew[0] = particle2frame1[1,0]
    Vnew[1] = particle2frame1[1,1]
    Vnew[2] = particle2frame1[1,2]
    Vnew[3] = particle2frame1[1,3]
    Jnew = np.zeros(4)
    Jnew[0] = particle2frame1[2,0]
    Jnew[1] = particle2frame1[2,1]
    Jnew[2] = particle2frame1[2,2]
    Jnew[3] = particle2frame1[2,3]

    particle1frame2 = FrameHop(particle2, particle1)
    xnew = np.zeros(4)
    xnew[0] = particle1frame2[0,0]
    xnew[1] = particle1frame2[0,1]
    xnew[2] = particle1frame2[0,2]
    xnew[3] = particle1frame2[0,3]
    vnew = np.zeros(4)
    vnew[0] = particle1frame2[1,0]
    vnew[1] = particle1frame2[1,1]
    vnew[2] = particle1frame2[1,2]
    vnew[3] = particle1frame2[1,3]
    jnew = np.zeros(4)
    jnew[0] = particle1frame2[2,0]
    jnew[1] = particle1frame2[2,1]
    jnew[2] = particle1frame2[2,2]
    jnew[3] = particle1frame2[2,3]

    gamma_factor = 1/np.sqrt(1-((Vnew[1] ** 2 + Vnew[2] ** 2 + Vnew[3] ** 2)/(Vnew[0]**2)))
    Gamma_Factor = 1/np.sqrt(1-((vnew[1] ** 2 + vnew[2] ** 2 + vnew[3] ** 2)/(vnew[0]**2)))
    n = (gamma_factor+Gamma_Factor) * np.sqrt(Xnew[1] ** 2 + Xnew[2] ** 2 + Xnew[3] ** 2)
    N = (gamma_factor+Gamma_Factor) * np.sqrt(xnew[1] ** 2 + xnew[2] ** 2 + xnew[3] ** 2)

    particle2frame1New = RetardedTime(m, a, Xnew, Vnew, Jnew, math.floor(n))
    particle1frame2New = RetardedTime(M, A, xnew, vnew, jnew, math.floor(N))

    #particle1New = FrameHop(particle1frame2New, particle2)
    #particle2New = FrameHop(particle2frame1New, particle1)
    particle1New = FrameHopInverse(particle2,particle1frame2New)
    particle2New = FrameHopInverse(particle1,particle2frame1New)

    particle1 = particle1New
    particle2 = particle2New

    data0 = np.array((particle2,particle1))
    return data0

def OffsetPair(data0):
    X = [data0[0,0],data0[0,1],data0[0,2],data0[0,3]]
    V = [data0[0,4],data0[0,5],data0[0,6],data0[0,7]]
    J = [data0[0,8],data0[0,9],data0[0,10],data0[0,11]]
    x = [data0[1,0],data0[1,1],data0[1,2],data0[1,3]]
    v = [data0[1,4],data0[1,5],data0[1,6],data0[1,7]]
    j = [data0[1,8],data0[1,9],data0[1,10],data0[1,11]]
    data1 = LinearPair(X, V, J, x, v, j)
    data2 = OrbitPair(X, V, J, x, v, j)
    data0 = data2-data1
    return data0

def Nbody(dataN,N):
    # print(dataN[1,:])
    data0 = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
    dataN0 = dataN * 0
    for n1 in range(0, N - 1):
        for n2 in range(n1 + 1, N):
            # print([n1,n2])
            data0[0, :] = dataN[n1, :]
            data0[1, :] = dataN[n2, :]
            data1 = OffsetPair(data0)
            dataN0[n1, :] += data1[0, :]
            dataN0[n2, :] += data1[1, :]
    dataN1 = dataN * 0
    for N in range(0, N):
        # print(N)
        data0 = dataN[N, :]
        X = [data0[0], data0[1], data0[2], data0[3]]
        V = [data0[4], data0[5], data0[6], data0[7]]
        J = [data0[8], data0[9], data0[10], data0[11]]
        data0 = LinearPair(X, V, J, [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
        data0 = data0[0, :]
        dataN1[N] = data0
    dataN = dataN1 + dataN0
    return dataN

def KerrOrbitPair(data0):
    X = [data0[0,0],data0[0,1],data0[0,2],data0[0,3]]
    V = [data0[0,4],data0[0,5],data0[0,6],data0[0,7]]
    J = [data0[0,8],data0[0,9],data0[0,10],data0[0,11]]
    x = [data0[1,0],data0[1,1],data0[1,2],data0[1,3]]
    v = [data0[1,4],data0[1,5],data0[1,6],data0[1,7]]
    j = [data0[1,8],data0[1,9],data0[1,10],data0[1,11]]
    data0 = OrbitPair(X, V, J, x, v, j)
    return data0

