
N = 12
n = 0
for iter0 in range(0,N+1):
    t0 = iter0 * 5
    for iter1 in range(0,iter0+1):
        r0 = t0 - 2 * iter1
        for iter2 in range(0,(iter0-iter1)+1):
            m0 = r0 * 3 - t0 * 2 - iter2 * 10
            n = n + 1
            #print([t0,r0,m0])
            c0 = t0 - r0
            G0 = 3 * r0 - m0 - 2*t0
            h0 = 2 * r0 + m0 - t0
            print([t0,r0,m0],[int(c0/2),int(G0/10),int(h0/10)])
print(n)
