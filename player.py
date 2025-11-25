import numpy as np
from blackhole import LorentzTransform
import pygame
from blackhole import Boost, Rotate

def KeyInput():
    keys = pygame.key.get_pressed()
    Ax = 0
    if keys[pygame.K_d]:
        Ax = 1
    elif keys[pygame.K_a]:
        Ax = -1
    Ay = 0
    if keys[pygame.K_w]:
        Ay = 1
    elif keys[pygame.K_s]:
        Ay = -1
    Az = 0
    if keys[pygame.K_r]:
        Az = 1
    elif keys[pygame.K_f]:
        Az = -1
    Tx = 0
    if keys[pygame.K_x]:
        Tx = 1
    elif keys[pygame.K_v]:
        Tx = -1
    Ty = 0
    if keys[pygame.K_q]:
        Ty = 1
    elif keys[pygame.K_e]:
        Ty = -1
    Tz = 0
    if keys[pygame.K_c]:
        Tz = 1
    elif keys[pygame.K_z]:
        Tz = -1
    dV6 = np.array((Ax,Ay,Az,Tx,Ty,Tz))/60
    return dV6

def StepMove(data,key0):
    dVx = key0[0]
    dVy = key0[1]
    dVz = key0[2]
    dJx = key0[3]
    dJy = key0[4]
    dJz = key0[5]
    t = data[0]
    x = data[1]
    y = data[2]
    z = data[3]
    dt = data[4]
    dx = data[5]
    dy = data[6]
    dz = data[7]
    Xx = data[8]
    Xy = data[9]
    Xz = data[10]
    Yx = data[11]
    Yy = data[12]
    Yz = data[13]
    Zx = data[14]
    Zy = data[15]
    Zz = data[16]
    Jx = data[17]
    Jy = data[18]
    Jz = data[19]
    dVx0 = dVx * Xx + dVy * Yx + dVz * Zx
    dVy0 = dVx * Xy + dVy * Yy + dVz * Zy
    dVz0 = dVx * Xz + dVy * Yz + dVz * Zz
    X = np.array((t,x,y,z))
    V = np.array((dt,dx,dy,dz))
    N = np.array([[1, 0, 0, 0],
                  [0 , Xx, Xy, Xz],
                  [0 , Yx, Yy, Yz],
                  [0 , Zx, Zy, Zz]])
    X = X + V
    Mboost = LorentzTransform(dVx0,dVy0,dVz0,0,0,0)
    V0 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            V0[u] += V[v]*Mboost[u,v]
    V = V0
    dJx0 = dJx * Xx + dJy * Yx + dJz * Zx
    dJy0 = dJx * Xy + dJy * Yy + dJz * Zy
    dJz0 = dJx * Xz + dJy * Yz + dJz * Zz
    Jx += dJx0
    Jy += dJy0
    Jz += dJz0
    Mrotate = LorentzTransform(0,0,0,Jx,Jy,Jz)
    N0 = np.zeros((4,4))
    for u in range(0, 4):
        for v in range(0, 4):
            for o in range(0,4):
                N0[u,o] += N[v,o] * Mrotate[u, v]
    N = N0
    data0 = np.array((X[0],X[1],X[2],X[3],V[0],V[1],V[2],V[3],N[1,1],N[1,2],N[1,3],N[2,1],N[2,2],N[2,3],N[3,1],N[3,2],N[3,3],Jx,Jy,Jz))
    return data0

def Camera(data, particle):
    t = data[0]
    x = data[1]
    y = data[2]
    z = data[3]
    dt = data[4]
    dx = data[5]
    dy = data[6]
    dz = data[7]
    Xx = data[8]
    Xy = data[9]
    Xz = data[10]
    Yx = data[11]
    Yy = data[12]
    Yz = data[13]
    Zx = data[14]
    Zy = data[15]
    Zz = data[16]
    Jx = data[17]
    Jy = data[18]
    Jz = data[19]
    t0 = particle[0]
    x0 = particle[1]
    y0 = particle[2]
    z0 = particle[3]
    t0 = t0-t
    x0 = x0-x
    y0 = y0-y
    z0 = z0-z
    Mboost = Boost(np.array((dt, dx, dy, dz)))
    X = np.array((t0, x0, y0, z0))
    X1 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += Mboost[u,v]*X[v]
    X = X1
    x1 = X[1]
    y1 = X[2]
    z1 = X[3]
    R2 = x1*x1 + y1*y1 + z1*z1
    X[0] = X[0] + np.sqrt(R2)
    MrotateZ1 = Rotate(np.array((0,Zx,Zy,Zz)))
    Zorient0 = [0,Xz,Yz,Zz]
    Xorient0 = [0, Xx, Yx, Zx]
    Yorient0 = [0, Xy, Yy, Zy]
    X1 = np.zeros(4)
    Xorient1 = np.zeros(4)
    Yorient1 = np.zeros(4)
    Zorient1 = np.zeros(4)
    for u in range(0,4):
        for v in range(0,4):
            X1[u] += MrotateZ1[u,v]*X[v]
            Xorient1[u] += MrotateZ1[u, v] * Xorient0[v]
            Yorient1[u] += MrotateZ1[u, v] * Yorient0[v]
            Zorient1[u] += MrotateZ1[u, v] * Zorient0[v]
    Xorient0 = Xorient1
    Yorient0 = Yorient1
    X = X1
    #print(Zorient1)
    Mrotatezx0 = Rotate(np.array((0,1,0,0)))
    Xorient1 = np.zeros(4)
    Yorient1 = np.zeros(4)
    for u in range(0, 4):
        for v in range(0, 4):
            X1[u] += Mrotatezx0[u, v] * X[v]
            Xorient1[u] += Mrotatezx0[u, v] * Xorient0[v]
            Yorient1[u] += Mrotatezx0[u, v] * Yorient0[v]
    X = X1
    MrotateXz1 = Rotate(np.array((0,Xorient1[1],Xorient1[2],Xorient1[3])))
    Xorient0 = Xorient1
    Xorient1 = np.zeros(4)
    Yorient0 = Yorient1
    Yorient1 = np.zeros(4)
    for u in range(0, 4):
        for v in range(0, 4):
            X1[u] += MrotateXz1[u, v] * X[v]
            Xorient1[u] += MrotateXz1[u, v] * Xorient0[v]
            Yorient1[u] += MrotateXz1[u, v] * Yorient0[v]
    MrotateZx1 = np.linalg.inv(Mrotatezx0)
    X = X1
    X1 = np.zeros(4)
    Xorient0 = Xorient1
    Xorient1 = np.zeros(4)
    Yorient0 = Yorient1
    Yorient1 = np.zeros(4)
    for u in range(0, 4):
        for v in range(0, 4):
            X1[u] += MrotateZx1[u,v]*X[v]
            Xorient1 += MrotateZx1[u,v]*Xorient0[v]
            Yorient1 += MrotateZx1[u,v]*Yorient0[v]
    #print(Xorient1)
    #print(Yorient1)
    X = X1
    t1 = X[0]
    x1 = X[1]
    y1 = X[2]
    z1 = X[3]

    m = particle[4]
    if (-12)<t1<0:
        xyz = [m,x1,y1,z1]
    else:
        xyz = [0,0,-1,0]
    #xyz = [m,0,1,0]
    return xyz