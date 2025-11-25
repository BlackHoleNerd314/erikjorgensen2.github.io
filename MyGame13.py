import numpy as np
import pygame
from pygame.locals import QUIT
import blackhole
import player

N = 800
n = 100

R = np.linspace(-n / 2, n / 2, n, endpoint=False)
X, Y, Z = np.meshgrid(R, R, R, indexing='ij')
r = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
Rscale0 = 1
Tscale0 = 1
Mscale0 = 1
pygame.init()
vector_norm = 1
screen = pygame.display.set_mode((N, N))
clock = pygame.time.Clock()
canvas = pygame.Surface((N, N))
position = np.array([0,50,0,0])
vector = np.array([1,0,0.1,0])
matrix = np.eye(4)
running = True
while running:
    for ev in pygame.event.get():
        if ev.type == QUIT:
            running = False
    scale0 = player.ScaleKeys()
    TimeScale = scale0[0] / 60
    SpaceScale = scale0[1] / 60
    MassScale = scale0[2] / 60

    dV6 = player.KeyInput()
    Vscale0 = Rscale0/Tscale0
    Matrix0 = blackhole.LorentzTransform(dV6[0]*Vscale0,dV6[1]*Vscale0,dV6[2]*Vscale0,dV6[3],dV6[4],dV6[5])
    matrix = np.matmul(matrix,Matrix0)
    vector = np.matmul(vector,Matrix0)
    Mscale0 = Mscale0 * np.exp(MassScale)
    #iter0 = iter1
    #data = blackhole.Geodesic(Mscale0 * iter0, 0, position * iter0, vector, [0, 0, 0, 0], np.int_(iter0*vector_norm)+5) / iter0
    #vector0 = (data[1, 0:4] - data[0, 0:4]) * iter0
    #vector_norm = vector[0]/vector0[0]
    #position = data[np.int_(iter0), 0:4]
    #vector31 = (data[np.int_(iter0)+1, 0:4] - data[np.int_(iter0)-1, 0:4]) * iter0/2
    #vector = vector_norm*vector31
    iter0 = np.ceil(Tscale0)
    data = blackhole.Geodesic(Mscale0/2,(1/2)/(Mscale0),position,vector,[0,0,0,0],np.int_(iter0)*3)
    position = data[np.int_(iter0), 0:4]
    vector = (data[np.int_(iter0)*2, 0:4] - data[0, 0:4])
    #print(position)



    vector = vector * np.exp(TimeScale)
    Tscale0 = Tscale0 * np.exp(TimeScale)
    Rscale0 = Rscale0 * np.exp(SpaceScale)
    X = X * np.exp(SpaceScale)
    Y = Y * np.exp(SpaceScale)
    Z = Z * np.exp(SpaceScale)
    r = r * np.exp(SpaceScale)

    def field():
        # Define 3D grid in player (local) coordinates
        flat_X = X.flatten()
        flat_Y = Y.flatten()
        flat_Z = Z.flatten()
        flat_r = r.flatten()

        # Retarded time in player's local frame
        flat_T_local = -flat_r  # observer sees emission from radius r ago
        events_local = np.vstack((flat_T_local, flat_X, flat_Y, flat_Z))  # Shape: (4, n^3)

        # Transform from player's local frame to world coordinates
        events_world = matrix @ events_local  # Apply Lorentz matrix
        events_world += position[:, None]  # Translate to global position

        # Extract world coordinates for field sampling
        t0 = events_world[0].reshape((n, n, n))
        x0 = events_world[1].reshape((n, n, n))
        y0 = events_world[2].reshape((n, n, n))
        z0 = events_world[3].reshape((n, n, n))

        # Sample world field — defined as static Gaussian in world space
        Psi = np.exp(-(x0 ** 2 + y0 ** 2 + z0 ** 2) / (2 * Mscale0**2))
        return Psi



    Psi = field()

    density = np.abs(Psi)**2

    density0 = np.sum(density,axis=2)

    # --- Create an 8‑bit grayscale RGB view --
    img = (density0*32).astype(np.uint8)

    # Convert to RGB grayscale (each pixel: (val, val, val))
    gray_rgb = np.stack([img] * 3, axis=-1)

    # Now make surface properly
    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (N, N))
    screen.blit(surf, (0, 0))
    pygame.display.flip()

    clock.tick(30)

pygame.quit()
