import numpy as np
import Update

def Trajectory(M,t,x,y,z,vx,vy,vz):
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    G0 = v*v*r/M
    h0 = v*r*M
    c0 = 1/v
    alpha1 = 137.036
    gravity0 = G0 < alpha1
    quantum0 = h0 < alpha1
    relativity0 = c0 < alpha1
    o1 = (not gravity0) and (not quantum0) and (not relativity0)
    C1 = (not gravity0) and (not quantum0) and (relativity0)
    G1 = (gravity0) and (not quantum0) and (not relativity0)
    H1 = (not gravity0) and (quantum0) and (not relativity0)
    CG1 = (gravity0) and (not quantum0) and (relativity0)
    CH1 = (not gravity0) and (quantum0) and (relativity0)
    GH1 = (gravity0) and (quantum0) and (not relativity0)
    CGH1 = (gravity0) and (quantum0) and (relativity0)
    if o1:
        return Update.Phys1(M,t,x,y,z,vx,vy,vz)
    if C1:
        return Update.PhysC(M, t, x, y, z, vx, vy, vz)
    if G1:
        return Update.PhysG(M,t,x,y,z,vx,vy,vz)
    if CG1:
        return Update.PhysCG(M,t,x,y,z,vx,vy,vz)
    if H1:
        return Update.PhysH(M,t,x,y,z,vx,vy,vz)
    if CH1:
        return Update.PhysCH(M,t,x,y,z,vx,vy,vz)
    if GH1:
        return Update.PhysGH(M,t,x,y,z,vx,vy,vz)
    if CGH1:
        return Update.PhysCGH(M,t,x,y,z,vx,vy,vz)


