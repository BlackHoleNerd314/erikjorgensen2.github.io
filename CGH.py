import numpy as np


tU = (4.35456 * 10**17)
c = 2.99792458 * 10**8
G = 6.67408 / 10**11
h = 6.62607015 / 10**34
kB = 1.380649 / 10**23

t0 = np.sqrt(h*G/c**5)
r0 = np.sqrt(h*G/c**3)
m0 = np.sqrt(h*c/G)

N = (tU/t0)

r1 = r0 * N
r2 = r0 * N**(3/5)

m1a = m0 * N
m1b = m0 / N
m2 = m0 / N**(1/5)

u1 = np.array([t0,r0,m0])
u2 = np.array([tU,r1,m1a])
u3 = np.array([tU,r1,m1b])
u4 = np.array([tU,r2,m2])
n1 = 1
n2 = 0
n3 = 0
n4 = 0
n0 = n1 + n2 + n3 + n4
u0 = u1**(n1/n0) * u2**(n2/n0) * u3**(n3/n0) * u4**(n4/n0)
print(u0,n0)
