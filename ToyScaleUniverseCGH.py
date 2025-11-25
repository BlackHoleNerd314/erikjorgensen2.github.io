import numpy as np
pi = np.pi
c = 2.99792458 * 10**8
G = 6.67408 / 10**11
h = 6.62607015 / 10**34
kB = 1.380649 / 10**23
C0 = 1/c
G0 = 8*pi*G
H0 = h/(2*pi)
K0 = kB
T0 = (2**200)*np.sqrt(H0*G0*(C0**5))
value0 = np.array([C0,G0,H0,K0,T0])
value1 = value0 ** 0.125
print(value1)


