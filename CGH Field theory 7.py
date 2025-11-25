import numpy as np
import pygame
from pygame.locals import *
screen_size = 800
x0 = 0
y0 = 0
dr = 100
dp = 0.01
G0 = 1
Lambda0 = -0
N = screen_size
A = 100
Mpsi = 1/screen_size**(1/3)
Mphi = 0/screen_size**(1/2)
R = np.linspace(-N/2, N/2, N, endpoint=False)
P = np.linspace(-1/2, 1/2, N, endpoint=False)
X, Y= np.meshgrid(R, R, indexing='ij')
Px, Py = np.meshgrid(P,P,indexing='ij')
noise = np.random.normal(loc=0.0, scale=1.0, size=X.shape)
Psihat = np.exp(-(Px**2 + Py**2)/(2*dp**2)) * noise
Psi = A * np.real(np.fft.ifft2(Psihat) * np.exp(-(X**2 + Y**2)/(2*dr**2)))
Phi = 0*Psi
PhiT = 0*Phi
PsiT = 0*Psi
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
    # Klein Gordon interaction loop
    density = np.abs(Psi) ** 2
    # Space Derivatives of Matter field
    dPsi_dx = (np.roll(Psi, -1, axis=0) - np.roll(Psi, 1, axis=0)) / 2.0
    dPsi_dy = (np.roll(Psi, -1, axis=1) - np.roll(Psi, 1, axis=1)) / 2.0
    d2Psi_dx2 = (np.roll(dPsi_dx, -1, axis=0) - np.roll(dPsi_dx, 1, axis=0)) / 2.0
    d2Psi_dy2 = (np.roll(dPsi_dy, -1, axis=1) - np.roll(dPsi_dy, 1, axis=1)) / 2.0
    # Matter Wave Loop(Runga Kutta Integration)
    LapPsi = (d2Psi_dx2 + d2Psi_dy2)
    PsiT = PsiT + LapPsi - Mpsi ** 2 * Psi + (Phi + Lambda0*density) * Psi
    Psi = Psi + PsiT
    # Space Derivatives of Force field
    dPhi_dx = (np.roll(Phi, -1, axis=0) - np.roll(Phi, 1, axis=0)) / 2.0
    dPhi_dy = (np.roll(Phi, -1, axis=1) - np.roll(Phi, 1, axis=1)) / 2.0
    d2Phi_dx2 = (np.roll(dPhi_dx, -1, axis=0) - np.roll(dPhi_dx, 1, axis=0)) / 2.0
    d2Phi_dy2 = (np.roll(dPhi_dy, -1, axis=1) - np.roll(dPhi_dy, 1, axis=1)) / 2.0
    # Force Wave Loop(Runga Kutta Integration)
    LapPhi = (d2Phi_dx2 + d2Phi_dy2)
    PhiT = PhiT + LapPhi - Mphi ** 2 * Phi + density * G0
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
    clock.tick(12)
pygame.quit()

