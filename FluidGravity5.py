import cupy as cp
import numpy as np
import pygame
from pygame.locals import *

# -----------------------------
# Parameters
# -----------------------------
N = 496  # number of particles
box_size = 100.0
cell_size = 3.0  # cutoff distance for Lennard-Jones
num_cells = int(box_size / cell_size)
dt = 0.001
steps = 100000

epsilon = 1.0
sigma = 1.0
G = 0.0001  # gravitational constant

# -----------------------------
# Particle initialization (2D uniform grid using linspace)
# -----------------------------
points_per_dim = int(cp.sqrt(N))
lin = cp.linspace(0, box_size, points_per_dim, endpoint=False)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N]
velocities = cp.zeros_like(positions)


# -----------------------------
# Helper: cell indexing
# -----------------------------
def get_cell_indices(pos):
    return cp.floor(pos / cell_size).astype(cp.int32)


def build_cell_list(positions):
    cell_list = {}
    indices = get_cell_indices(positions)
    for i, idx in enumerate(indices):
        key = tuple(idx.tolist())
        if key not in cell_list:
            cell_list[key] = []
        cell_list[key].append(i)
    return cell_list


# -----------------------------
# Lennard-Jones force with short-range velocity damping
# -----------------------------
def compute_lj_forces(positions, velocities):
    forces = cp.zeros_like(positions)
    cell_list = build_cell_list(positions)
    damping_coeff = 0.01  # small damping

    for cell, particles in cell_list.items():
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                neighbor = (cell[0] + dx, cell[1] + dy)
                if neighbor in cell_list:
                    neighbors.extend(cell_list[neighbor])

        for i in particles:
            for j in neighbors:
                if i < j:
                    rij = positions[i] - positions[j]
                    rij = rij - box_size * cp.rint(rij / box_size)  # periodic BC
                    r2 = cp.dot(rij, rij)
                    if r2 < cell_size ** 2 and r2 > 1e-12:
                        inv_r6 = (sigma ** 2 / r2) ** 3
                        f_scalar = 24 * epsilon * inv_r6 * (2 * inv_r6 - 1) / r2
                        fij = f_scalar * rij
                        forces[i] += fij
                        forces[j] -= fij

                        # Relative velocity damping
                        vij = velocities[i] - velocities[j]
                        damping_force = -damping_coeff * vij
                        forces[i] += damping_force
                        forces[j] -= damping_force
    return forces


# -----------------------------
# Gravitational force via Fourier pseudospectral Poisson solver
# -----------------------------
def compute_gravity_forces(positions):
    grid_size = num_cells
    rho = cp.zeros((grid_size, grid_size), dtype=cp.float32)

    # deposit mass to grid using nearest grid point
    idx = cp.floor(positions / box_size * grid_size).astype(cp.int32) % grid_size
    for i in range(len(positions)):
        rho[idx[i, 0], idx[i, 1]] += 1.0

    # Fourier transform
    rho_k = cp.fft.fft2(rho)
    kx = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0  # avoid division by zero

    phi_k = -4 * cp.pi * G * rho_k / K2
    phi = cp.fft.ifft2(phi_k).real

    # Compute gradient for forces
    fx_k = (1j * KX) * phi_k
    fy_k = (1j * KY) * phi_k
    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real

    forces = cp.zeros_like(positions)
    for i in range(len(positions)):
        gx, gy = idx[i]
        forces[i, 0] = Fx[gx, gy]
        forces[i, 1] = Fy[gx, gy]
    return forces


# -----------------------------
# Unified force computation
# -----------------------------
def compute_forces(positions, velocities):
    forces_lj = compute_lj_forces(positions, velocities)
    forces_g = compute_gravity_forces(positions)
    return forces_lj + forces_g


# -----------------------------
# Velocity-Verlet integration
# -----------------------------
def integrate(positions, velocities, forces):
    positions += velocities * dt + 0.5 * forces * dt ** 2
    positions %= box_size
    new_forces = compute_forces(positions, velocities)
    velocities += 0.5 * (forces + new_forces) * dt
    return positions, velocities, new_forces


# -----------------------------
# Pygame visualization
# -----------------------------
pygame.init()
screen = pygame.display.set_mode((600, 600))
clock = pygame.time.Clock()

forces = compute_forces(positions, velocities)

for step in range(steps):
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            exit()

    positions, velocities, forces = integrate(positions, velocities, forces)

    screen.fill((0, 0, 0))
    pos_np = cp.asnumpy(positions / box_size * 600).astype(int)
    for x, y in pos_np:
        pygame.draw.circle(screen, (255, 255, 255), (x, y), 2)
    pygame.display.flip()
    clock.tick(60)
