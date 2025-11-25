import numpy as np
from itertools import permutations

def FourierSeries(N):
    X = np.linspace(0,1-1/N,N)
    P = np.linspace(0,1-1/N,N)
    x , p = np.meshgrid(X, P, indexing='ij')
    phase0 = np.exp(-2j*x*p*np.pi*N)
    return phase0

def count_inversions(p):
    # Simple insertion sort-style inversion counter
    inv = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv += 1
    return 1 if inv % 2 == 0 else -1

def generate_permutations_with_parity(N):
    perm_list = list(permutations(range(1, N + 1)))
    perms = np.array(perm_list, dtype=np.int32)
    signs = np.array([count_inversions(p) for p in perms], dtype=np.int8)
    return np.column_stack((perms, signs))

def Factorial0(N):
    n0 = 1
    for n in range(1,N+1):
        n0 = n0 * n
    return n0

def Entangle(N):
    basis = FourierSeries(N)
    perm0 = generate_permutations_with_parity(N)
    det = 0
    prm = 0
    for n in range(0,Factorial0(N)):
        list0 = perm0[n,0:N]
        sign0 = perm0[n,N]
        det0 = sign0
        prm0 = 1
        for n1 in range(0,N):
            n2 = list0[n1] - 1
            det0 = det0 * basis[n1,n2]
            prm0 = prm0 * basis[n1,n2]
        det = det + det0
        prm = prm + prm0
    return det,prm

