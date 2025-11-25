import numpy as np
import pygame
from pygame.locals import *

# -----------------------------
# Simulation parameters
# -----------------------------
M = 0.25
G = 0.001
m = 0
x0, y0, z0 = 5, -2, 3
dr = 4
px0, py0, pz0 = 0.15, -0.10, -0.06

# Grid size
screen_size = N = int(1 / M**3)  # 64
R = np.linspace(-N/2, N/2, N, endpoint=False)
X, Y, Z = np.meshgrid(R, R, R, indexing='ij')

# Initialize fields
Psi = np.exp(-((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)/(2*dr**2)) * np.exp(-1j*(px0*X+py0*Y+pz0*Z))
PsiT = np.zeros_like(Psi)
Phi = np.zeros_like(Psi, dtype=np.float64)
PhiT = np.zeros_like(Psi, dtype=np.float64)

# -----------------------------
# Pygame setup
# -----------------------------
pygame.init()
screen = pygame.display.set_mode((screen_size, screen_size))
clock = pygame.time.Clock()

# -----------------------------
# Helper function
# -----------------------------
def normalize_density(density):
    density = np.real(density)
    # Proper normalization to [0, 255]
    density_norm = (density - density.min()) / (density.max() - density.min()) * 255
    return density_norm.astype(np.uint8)

# -----------------------------
# Main simulation loop
# -----------------------------
running = True
dt = 0.1  # time step

while running:
    for event in pygame.event.get():
        if event.type == QUIT:
            running = False

    # -----------------------------
    # Compute Laplacians
    # -----------------------------
    LapPsi = (
        np.roll(Psi, -1, axis=0) - 2*Psi + np.roll(Psi, 1, axis=0) +
        np.roll(Psi, -1, axis=1) - 2*Psi + np.roll(Psi, 1, axis=1) +
        np.roll(Psi, -1, axis=2) - 2*Psi + np.roll(Psi, 1, axis=2)
    )

    LapPhi = (
        np.roll(Phi, -1, axis=0) - 2*Phi + np.roll(Phi, 1, axis=0) +
        np.roll(Phi, -1, axis=1) - 2*Phi + np.roll(Phi, 1, axis=1) +
        np.roll(Phi, -1, axis=2) - 2*Phi + np.roll(Phi, 1, axis=2)
    )

    # -----------------------------
    # Update fields (Euler step)
    # -----------------------------
    PsiT += dt * (LapPsi - M**2 * Psi + G * Psi * Phi - Psi * np.abs(Psi)**2)
    Psi += dt * PsiT

    PhiT += dt * (LapPhi - m**2 * Phi + np.abs(Psi)**2)
    Phi += dt * PhiT

    #---------------
    #Remove constant component to the massless fields
    #---------------

    Phi0 = (np.sum(Phi)/N**3)*np.ones([N,N,N])
    Phi -= Phi0

    # -----------------------------
    # Visualization (sum over z-axis)
    # -----------------------------
    density2d = np.sum(np.abs(Psi)**2, axis=2)
    img1 = normalize_density(density2d)
    density2d = np.sum(np.abs(Phi)**2, axis=2)
    img0 = normalize_density(density2d)
    gray_rgb = np.stack([img1,img0,0*img0], axis=-1)

    # Convert to Pygame surface
    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (screen_size, screen_size))
    screen.blit(surf, (0, 0))
    pygame.display.flip()
    clock.tick(60)

pygame.quit()
