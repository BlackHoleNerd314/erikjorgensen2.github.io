import numpy as np
from scipy.linalg import expm
import BlackHoleframe

i = 1j

def J(s,a0,a):
    Jx = 0
    Jy = 0
    Jz = 0
    term0 = 0
    term1 = 0
    if a0==a+1:
        term0 = (s-a)*(s+a+1)
    if a0==a-1:
        term1 = (s+a)*(s-a+1)
    Jx = Jx + np.sqrt(term0)
    Jx = Jx + np.sqrt(term1)
    Jy = Jy + np.sqrt(term0)
    Jy = Jy - np.sqrt(term1)
    Jx = Jx / 2
    Jy = Jy / (2*i)
    if a0==a:
        Jz = a
    return Jx,Jy,Jz

def SpinMatrix(s):
    Jx = np.complex128(np.zeros((s, s)))
    Jy = np.complex128(np.zeros((s, s)))
    Jz = np.complex128(np.zeros((s, s)))
    for s0 in range(0,s):
        for s1 in range(0,s):
            j0 = (s-1)/2
            x0,y0,z0 = J(j0,s0-j0,s1-j0)
            Jx[s0, s1] = np.complex128(x0)
            Jy[s0, s1] = np.complex128(y0)
            Jz[s0, s1] = np.complex128(z0)
    return Jx,Jy,Jz

def Psi(m,a,t,x,y,z,Kx,Ky,Kz,Jx,Jy,Jz):
    s = 2*(a*m) + 1
    Jx0,Jy0,Jz0 = SpinMatrix(s)
    Mx = (Kx+i*Jx)*Jx0
    My = (Ky+i*Jy)*Jy0
    Mz = (Kz+i*Jz)*Jz0
    M = Mx+My+Mz
    M0 = expm(M)
    return M0











