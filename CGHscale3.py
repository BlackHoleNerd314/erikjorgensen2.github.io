import numpy as np
n = 0
N = 60
for t in range(0,N+1):
    r0 = int(np.floor(3*t/5))
    for r in range(r0,t+1):
        m1 = int(np.floor(r*3-t*2))
        m2 = int(np.floor(t*1-r*2))
        for m in range(m2,m1+1):
            C = t-r
            G = r*3-t*2-m
            H = r*2-t+m
            print([t-43,r-35,m-8],[C,G,H])
            n = n + 1
print(n)