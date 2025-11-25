import numpy as np
import pygame
from pygame.locals import *
M = 0.25
screen_size = np.int32(1/M**3)
N = screen_size
x0 = 5
y0 = -2
z0 = 3
dr = 4
px0 = 0.15
py0 = -0.10
pz0 = -0.06
G = 0.015625
R = np.linspace(-N/2, N/2, N, endpoint=False)
X, Y, Z= np.meshgrid(R, R, R, indexing='ij')
Psi = 0.1 * np.real(np.exp(-((X-x0)**2 + (Y-y0)**2 + (Z-z0)**2)/(2*dr**2))*np.exp(-1j*(px0*X+py0*Y+pz0*Z)))
PsiT = 0*Psi
Phi = 0*Psi
PhiT = 0*Psi
# Initialize Pygame
pygame.init()
screen = pygame.display.set_mode((screen_size, screen_size))
clock = pygame.time.Clock()
def normalize_density(density):
    """Normalize density for visualization"""
    density = np.real(density)
    return ((density+np.min(density)) / (np.max(density) + np.min(density))) * 255
running = True
while running:
    for event in pygame.event.get():
        if event.type == QUIT:
            running = False
    # Space Derivatives of Matter field
    d_dx = (np.roll(Psi, -1, axis=0) - np.roll(Psi, 1, axis=0)) / 2.0
    d_dy = (np.roll(Psi, -1, axis=1) - np.roll(Psi, 1, axis=1)) / 2.0
    d_dz = (np.roll(Psi, -1, axis=2) - np.roll(Psi, 1, axis=2)) / 2.0
    d2_dx2 = (np.roll(d_dx, -1, axis=0) - np.roll(d_dx, 1, axis=0)) / 2.0
    d2_dy2 = (np.roll(d_dy, -1, axis=1) - np.roll(d_dy, 1, axis=1)) / 2.0
    d2_dz2 = (np.roll(d_dz, -1, axis=2) - np.roll(d_dz, 1, axis=2)) / 2.0
    # Matter Wave Loop(Runga Kutta Integration)
    LapPsi = (d2_dx2 + d2_dy2 + d2_dz2)
    PsiT = PsiT + LapPsi - M**2*Psi + G*Psi*Phi
    Psi = Psi + PsiT
    # Space Derivatives of Force field
    d_dx = (np.roll(Phi, -1, axis=0) - np.roll(Phi, 1, axis=0)) / 2.0
    d_dy = (np.roll(Phi, -1, axis=1) - np.roll(Phi, 1, axis=1)) / 2.0
    d_dz = (np.roll(Phi, -1, axis=2) - np.roll(Phi, 1, axis=2)) / 2.0
    d2_dx2 = (np.roll(d_dx, -1, axis=0) - np.roll(d_dx, 1, axis=0)) / 2.0
    d2_dy2 = (np.roll(d_dy, -1, axis=1) - np.roll(d_dy, 1, axis=1)) / 2.0
    d2_dz2 = (np.roll(d_dz, -1, axis=2) - np.roll(d_dz, 1, axis=2)) / 2.0
    # Force Wave Loop(Runga Kutta Integration)
    LapPhi = (d2_dx2 + d2_dy2 + d2_dz2)
    PhiT = PhiT + LapPhi + Psi*np.conj(Psi)
    Phi = Phi + PhiT
    # Visualization
    density0 = np.sum(np.abs(Psi)**2, axis=2)
    img = normalize_density(density0).astype(np.uint8)
    # Convert to RGB grayscale (each pixel: (val, val, val))
    gray_rgb = np.stack([img] * 3, axis=-1)
    # Now make surface properly
    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (screen_size, screen_size))
    screen.blit(surf, (0, 0))
    pygame.display.flip()
    clock.tick(60)
pygame.quit()

