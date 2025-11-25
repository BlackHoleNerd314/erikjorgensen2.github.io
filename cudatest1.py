import cupy as cp
print(cp.__version__)
import numpy as np  # still needed for Pygame conversion
import pygame
from pygame.locals import *

# Parameters
M = 0.2
screen_size = int(1/M**3)
N = screen_size
x0, y0, z0 = 5, -2, 3
dr = 10
px0, py0, pz0 = 0.15, -0.10, -0.06
G = 0.001
contact0 = 0

# GPU arrays
R = cp.linspace(-N/2, N/2, N, endpoint=False)
X, Y, Z = cp.meshgrid(R, R, R, indexing='ij')
Psi = cp.real(cp.exp(-((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)/(2*dr**2)) *
              cp.exp(-1j*(px0*X+py0*Y+pz0*Z)))
PsiT = cp.zeros_like(Psi)
Phi = cp.zeros_like(Psi)
PhiT = cp.zeros_like(Psi)

# Pygame init
pygame.init()
screen = pygame.display.set_mode((screen_size, screen_size))
clock = pygame.time.Clock()

def normalize_density(density):
    """Normalize density for visualization"""
    density = cp.real(density)
    min_val = cp.min(density)
    max_val = cp.max(density)
    return ((density - min_val) / (max_val - min_val + 1e-9)) * 255

running = True
while running:
    for event in pygame.event.get():
        if event.type == QUIT:
            running = False

    # Space Derivatives of Matter field
    d_dx = (cp.roll(Psi, -1, axis=0) - cp.roll(Psi, 1, axis=0)) / 2.0
    d_dy = (cp.roll(Psi, -1, axis=1) - cp.roll(Psi, 1, axis=1)) / 2.0
    d_dz = (cp.roll(Psi, -1, axis=2) - cp.roll(Psi, 1, axis=2)) / 2.0

    d2_dx2 = (cp.roll(d_dx, -1, axis=0) - cp.roll(d_dx, 1, axis=0)) / 2.0
    d2_dy2 = (cp.roll(d_dy, -1, axis=1) - cp.roll(d_dy, 1, axis=1)) / 2.0
    d2_dz2 = (cp.roll(d_dz, -1, axis=2) - cp.roll(d_dz, 1, axis=2)) / 2.0

    LapPsi = d2_dx2 + d2_dy2 + d2_dz2
    PsiT = PsiT + LapPsi - M**2 * Psi + G * Psi * Phi - (Psi * cp.conj(Psi)) * contact0 * Psi
    Psi = Psi + PsiT

    # Space Derivatives of Force field
    d_dx = (cp.roll(Phi, -1, axis=0) - cp.roll(Phi, 1, axis=0)) / 2.0
    d_dy = (cp.roll(Phi, -1, axis=1) - cp.roll(Phi, 1, axis=1)) / 2.0
    d_dz = (cp.roll(Phi, -1, axis=2) - cp.roll(Phi, 1, axis=2)) / 2.0

    d2_dx2 = (cp.roll(d_dx, -1, axis=0) - cp.roll(d_dx, 1, axis=0)) / 2.0
    d2_dy2 = (cp.roll(d_dy, -1, axis=1) - cp.roll(d_dy, 1, axis=1)) / 2.0
    d2_dz2 = (cp.roll(d_dz, -1, axis=2) - cp.roll(d_dz, 1, axis=2)) / 2.0

    LapPhi = d2_dx2 + d2_dy2 + d2_dz2
    PhiT = PhiT + LapPhi + Psi * cp.conj(Psi)
    Phi = Phi + PhiT

    # Visualization (reduce to 2D and bring back to CPU for Pygame)
    density0 = cp.sum(cp.abs(Psi)**2, axis=2)
    img = normalize_density(density0).astype(cp.uint8)

    # Transfer to CPU (costly, so only 2D slice)
    img_cpu = cp.asnumpy(img)

    # Convert to RGB grayscale
    gray_rgb = np.stack([img_cpu] * 3, axis=-1)
    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (screen_size, screen_size))
    screen.blit(surf, (0, 0))

    pygame.display.flip()
    clock.tick(60)

pygame.quit()
