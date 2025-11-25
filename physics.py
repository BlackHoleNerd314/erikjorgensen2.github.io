import numpy as np

def Scale(m,kg,s,K):
    C = 100
    cycle = 1
    mol = 1
    c = 2.99792458*10**8
    G = 6.67408/10**11
    h = 6.62607015/10**34
    c0 = m/s
    G0 = (m*m*m)/(kg*s*s)
    h0 = (m*m*kg)/s
    kB = 1.380649/10**23
    kB0 = (m*m*kg)/(s*s*K)
    Ke = c**2/10**7
    Ke0 = G0*((kg/C)**2)
    Qe = 1.602176634/10**19
    Me = (1/10**3)*(5.485799090441/10**4)/(6.02214076*10**23)
    Qe0 = C
    Me0 = kg
    h0 = h0/(cycle*mol)
    kB0 = kB0/mol
    Qe0 = Qe0/mol
    Me0 = Me0/mol
    tUniverse = 4.35456*10**17
    tUniverse0 = s
    constants = [c/c0,G/G0,h/h0,kB/kB0,tUniverse/tUniverse0,Ke/Ke0,Qe/Qe0,Me/Me0]
    n1 = 1/(c/c0)
    n2 = 2*(G/G0)
    n3 = (h/h0)/(2*np.pi)
    n4 = 1/(tUniverse/tUniverse0)
    n5 = kB/kB0
    constants1 = [n1,n2,n3,n4,n5]
    return constants1

#constants = Scale(1609.344,0.001,3600,100,1)*np.ones((5,8))
constants1 = Scale(1,1,1,1)*np.ones((4,5))
unit0 = np.array([[1,0,-1,0,0],[3,-1,-2,0,0],[2,1,-1,0,0],[2,1,-2,-1,0],[0,0,1,0,0],[3,1,-2,0,-2],[0,0,0,0,1],[0,1,0,0,0]])
unit1 = np.array([[-1,0,1,0],[3,-1,-2,0],[2,1,-1,0],[0,0,-1,0],[2,1,-2,-1]])

print(np.transpose(constants1)**(1/unit1))