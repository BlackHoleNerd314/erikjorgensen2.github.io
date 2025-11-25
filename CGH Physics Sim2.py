import numpy as np
import pygame
from pygame.locals import *
screen_size = 500
x0 = 50
y0 = -20
dr = 35
px0 = 0.1
py0 = 0.3
G = 0.01
N = screen_size
M = 0.1
R = np.linspace(-N/2, N/2, N, endpoint=False)
P = np.linspace(-1/2, 1/2, N, endpoint=False)
X, Y= np.meshgrid(R, R, indexing='ij')
Px, Py= np.meshgrid(P, P, indexing='ij')
k_squared = Px**2 + Py**2
omega = np.sqrt(M**2 + k_squared)
K = np.sqrt(k_squared)
Psi = np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*dr**2))*np.exp(-1j*(px0*X+py0*Y))
Phi = 0*Psi
PhiT = 0*Phi
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
    # Klein Gordon matter Loop(Fourier Pseudo Spectral)
    PsiHat = np.fft.fft2(Psi)
    PsiHat = PsiHat * np.exp(1j * omega)
    Psi = np.fft.ifft2(PsiHat)
    # Klein Gordon interaction loop
    density = np.abs(Psi) ** 2
    Psi = Psi * np.exp(1j * Phi)
    Psi = Psi / np.sqrt(np.sum(np.sum(density)))
    # Space Derivatives of Force field
    d_dx = (np.roll(Phi, -1, axis=0) - np.roll(Phi, 1, axis=0)) / 2.0
    d_dy = (np.roll(Phi, -1, axis=1) - np.roll(Phi, 1, axis=1)) / 2.0
    d2_dx2 = (np.roll(d_dx, -1, axis=0) - np.roll(d_dx, 1, axis=0)) / 2.0
    d2_dy2 = (np.roll(d_dy, -1, axis=1) - np.roll(d_dy, 1, axis=1)) / 2.0
    # Force Wave Loop(Runga Kutta Integration)
    LapPhi = (d2_dx2 + d2_dy2)
    PhiT = PhiT + LapPhi - 0**2*Phi + density*G
    Phi = Phi + PhiT
    # Visualization
    density0 = density #np.sum(density, axis=2)
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

