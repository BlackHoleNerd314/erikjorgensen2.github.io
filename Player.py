import numpy as np
import pygame
import BlackHoleframe

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

def PlayerInit():
    PlayerPosition = np.array([0, 0, 0, 0])
    PlayerVelocity = np.array([1, 0, 0, 0])
    PlayerOrientation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    PlayerAngularMomentum = np.array([0, 0, 1])
    PlayerMass = 1
    return PlayerPosition, PlayerVelocity, PlayerOrientation, PlayerAngularMomentum, PlayerMass

def ScaleInit():
    TimeScale = 1
    SpaceScale = 1
    MassScale = 1
    return TimeScale, SpaceScale, MassScale

def PlayerScale(TimeScale, SpaceScale, MassScale):
    keys = pygame.key.get_pressed()
    if keys[pygame.K_t]:
        TimeScale *= 1.01
    elif keys[pygame.K_b]:
        TimeScale /= 1.01
    if keys[pygame.K_y]:
        SpaceScale *= 1.01
    elif keys[pygame.K_n]:
        SpaceScale /= 1.01
    if keys[pygame.K_u]:
        MassScale *= 1.01
    elif keys[pygame.K_m]:
        MassScale /= 1.01
    if keys[pygame.K_j]:
        MassScale = 1
    if keys[pygame.K_h]:
        SpaceScale = 1
    if keys[pygame.K_g]:
        TimeScale = 1
    return TimeScale, SpaceScale, MassScale

def PlayerUpdate(TimeScale, SpaceScale, MassScale, PlayerPosition, PlayerVelocity, PlayerOrientation, PlayerAngularMomentum, PlayerMass):
    TimeScale, SpaceScale, MassScale = PlayerScale(TimeScale, SpaceScale, MassScale)
    VelocityScale = SpaceScale / TimeScale
    dV6 = KeyInput()
    Jx = PlayerAngularMomentum[0]
    Jy = PlayerAngularMomentum[1]
    Jz = PlayerAngularMomentum[2]
    rotate0 = BlackHoleframe.LorentzTransform(0, 0, 0, Jx, Jy, Jz)
    PlayerOrientation = PlayerOrientation @ rotate0
    dV0 = np.array([dV6[0],dV6[1],dV6[2]])
    dv0 = dV0 @ PlayerOrientation
    Kx = dv0[0] * VelocityScale
    Ky = dv0[1] * VelocityScale
    Kz = dv0[2] * VelocityScale
    boostVel = BlackHoleframe.LorentzTransform(Kx,Ky,Kz,0,0,0)
    PlayerVelocity = PlayerVelocity @ boostVel
    PlayerPosition += PlayerVelocity * TimeScale
    dV1 = np.array([dV6[3], dV6[4], dV6[5]])
    dv1 = dV1 @ PlayerOrientation
    PlayerAngularMomentum[0] += dv1[0]
    PlayerAngularMomentum[1] += dv1[1]
    PlayerAngularMomentum[2] += dv1[2]
    return TimeScale, SpaceScale, MassScale, PlayerPosition, PlayerVelocity, PlayerOrientation, PlayerAngularMomentum, PlayerMass





























