import numpy as np
n = 0
N = 60
for t in range(0,N+1):
    r0 = int(np.floor(3*t/5))
    for r in range(r0,t+1):
        print([t-43,r-34])
        n = n + 1
print(n)