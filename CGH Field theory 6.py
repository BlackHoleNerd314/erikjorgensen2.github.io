import numpy as np
import pygame
from numba import njit, prange

# -----------------------------
# Simulation parameters
# -----------------------------
M = 0.25
G = 0.015625
x0, y0, z0 = 5, -2, 3
dr = 4
px0, py0, pz0 = 0.15, -0.10, -0.06
screen_size = N = int(1 / M ** 3)

R = np.linspace(-N / 2, N / 2, N, endpoint=False)
X, Y, Z = np.meshgrid(R, R, R, indexing='ij')

# Initialize fields
Psi = np.exp(-((X - x0) ** 2 + (Y - y0) ** 2 + (Z - z0) ** 2) / (2 * dr ** 2)) * np.exp(
    -1j * (px0 * X + py0 * Y + pz0 * Z))
PsiT = np.zeros_like(Psi)
Phi = np.zeros_like(Psi, dtype=np.float64)
PhiT = np.zeros_like(Psi, dtype=np.float64)
dt = 0.1


# -----------------------------
# CPU-optimized update kernel
# -----------------------------
@njit(parallel=True)
def update_fields(Psi, PsiT, Phi, PhiT, M, G, dt):
    N = Psi.shape[0]
    for i in prange(N):
        for j in range(N):
            for k in range(N):
                ip = (i + 1) % N
                im = (i - 1) % N
                jp = (j + 1) % N
                jm = (j - 1) % N
                kp = (k + 1) % N
                km = (k - 1) % N

                lap_Psi = (Psi[ip, j, k] + Psi[im, j, k] +
                           Psi[i, jp, k] + Psi[i, jm, k] +
                           Psi[i, j, kp] + Psi[i, j, km] -
                           6 * Psi[i, j, k])

                lap_Phi = (Phi[ip, j, k] + Phi[im, j, k] +
                           Phi[i, jp, k] + Phi[i, jm, k] +
                           Phi[i, j, kp] + Phi[i, j, km] -
                           6 * Phi[i, j, k])

                PsiT[i, j, k] += dt * (lap_Psi - M ** 2 * Psi[i, j, k] + G * Psi[i, j, k] * Phi[i, j, k])
                Psi[i, j, k] += dt * PsiT[i, j, k]

                PhiT[i, j, k] += dt * (lap_Phi + (Psi[i, j, k].real ** 2 + Psi[i, j, k].imag ** 2))
                Phi[i, j, k] += dt * PhiT[i, j, k]


# -----------------------------
# Pygame setup
# -----------------------------
pygame.init()
screen = pygame.display.set_mode((screen_size, screen_size))
clock = pygame.time.Clock()


def normalize_density(density):
    density_norm = (density - density.min()) / (density.max() - density.min()) * 255
    return density_norm.astype(np.uint8)


# -----------------------------
# Main loop
# -----------------------------
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    update_fields(Psi, PsiT, Phi, PhiT, M, G, dt)

    density2d = np.sum(np.abs(Psi) ** 2, axis=2)
    img = normalize_density(density2d)
    gray_rgb = np.stack([img] * 3, axis=-1)

    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (screen_size, screen_size))
    screen.blit(surf, (0, 0))
    pygame.display.flip()
    clock.tick(60)

pygame.quit()
