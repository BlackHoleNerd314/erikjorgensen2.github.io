import cupy as cp
import numpy as np
import pygame
from pygame.locals import *

# -----------------------------
# Parameters (same semantics as yours)
# -----------------------------
N = 512
box_size = 100.0
cell_size = 3.0   # LJ cutoff distance
num_cells = int(box_size / cell_size)
dt = 0.0005       # slightly smaller dt for stability; tune as needed
steps = 100000

epsilon = 1.0
sigma = 1.0
G = 0.0001

# -----------------------------
# Particle initialization (GPU arrays)
# -----------------------------
points_per_dim = int(cp.sqrt(N))
lin = cp.linspace(0, box_size, points_per_dim, endpoint=False, dtype=cp.float32)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N].astype(cp.float32)  # shape (N,2)
velocities = cp.zeros_like(positions)

# -----------------------------
# Vectorized Lennard-Jones + pairwise damping (GPU)
# -----------------------------
cutoff2 = float(cell_size ** 2)
damping_coeff = 0.01
eps_small = 1e-12

def compute_lj_forces(positions, velocities):
    """
    Fully vectorized LJ pairwise forces with periodic BC and short-range velocity damping.
    Returns array of shape (N,2) of forces on each particle.
    """
    # Pairwise displacement (N, N, 2)
    rij = positions[:, None, :] - positions[None, :, :]

    # Apply periodic BC (minimum image)
    rij = rij - box_size * cp.rint(rij / box_size)

    # Squared distances (N, N)
    r2 = cp.sum(rij * rij, axis=2)

    # Prevent self-interaction: set diagonal large
    r2 += cp.eye(N, dtype=r2.dtype) * 1e6

    # Neighbor mask for cutoff
    neighbor_mask = r2 < cutoff2  # shape (N, N)

    # Compute LJ scalar only where neighbor_mask is True
    # Avoid division by zero by using r2 + eps_small (diagonal already large)
    inv_r2 = 1.0 / (r2 + eps_small)
    inv_r6 = (sigma * sigma * inv_r2) ** 3  # (N,N)
    f_scalar = 24.0 * epsilon * inv_r6 * (2.0 * inv_r6 - 1.0) * inv_r2  # (N,N)

    # Zero out outside cutoff
    f_scalar = f_scalar * neighbor_mask.astype(f_scalar.dtype)  # (N,N)

    # Pairwise force vectors (N,N,2)
    fij = f_scalar[..., None] * rij

    # Sum over j to get total LJ force on i
    forces_lj = cp.sum(fij, axis=1)  # (N,2)

    # --- Damping: sum pairwise -damping_coeff*(v_i - v_j) over neighbors ---
    vij = velocities[:, None, :] - velocities[None, :, :]  # (N,N,2)
    neighbor_mask_3 = neighbor_mask[..., None]  # (N,N,1)
    damping_force = -damping_coeff * cp.sum(neighbor_mask_3 * vij, axis=1)  # (N,2)

    forces = forces_lj + damping_force
    return forces.astype(cp.float32)

# -----------------------------
# Vectorized gravitational PDE via FFT (GPU)
# -----------------------------
def compute_gravity_forces(positions):
    """
    Particle->grid mass deposit (nearest grid point), solve Poisson by FFT,
    compute gradient and sample forces at particle cell centers.
    """
    grid_size = num_cells
    rho = cp.zeros((grid_size, grid_size), dtype=cp.float32)

    # Nearest-grid-point deposit indices (grid coords)
    # Map positions in [0, box_size) to grid indices 0..grid_size-1
    idx = cp.floor(positions / box_size * grid_size).astype(cp.int32) % grid_size

    # GPU-side scatter add
    cp.add.at(rho, (idx[:, 0], idx[:, 1]), 1.0)

    # Fourier solve
    rho_k = cp.fft.fft2(rho)

    # wavevectors
    dk = box_size / grid_size
    kx = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')

    K2 = KX**2 + KY**2
    # avoid division by zero; set mean mode phi_k(0,0) = 0 (no net potential)
    K2[0, 0] = 1.0
    phi_k = -4.0 * cp.pi * G * rho_k / K2
    phi_k[0, 0] = 0.0

    # Gradient in k-space (i*k * phi_k)
    fx_k = (1j * KX) * phi_k
    fy_k = (1j * KY) * phi_k

    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real

    # Sample grid forces at particle grid locations
    forces = cp.empty_like(positions)
    gx = idx[:, 0] % grid_size
    gy = idx[:, 1] % grid_size
    forces[:, 0] = Fx[gx, gy]
    forces[:, 1] = Fy[gx, gy]

    return forces.astype(cp.float32)

# -----------------------------
# Unified force computation (GPU)
# -----------------------------
def compute_forces(positions, velocities):
    lj = compute_lj_forces(positions, velocities)
    grav = compute_gravity_forces(positions)
    return lj + grav

# -----------------------------
# Velocity-Verlet integration (GPU)
# -----------------------------
def integrate(positions, velocities, forces):
    # positions: (N,2), velocities: (N,2), forces: (N,2)
    positions = (positions + velocities * dt + 0.5 * forces * (dt ** 2)) % box_size
    new_forces = compute_forces(positions, velocities)
    velocities = velocities + 0.5 * (forces + new_forces) * dt
    return positions, velocities, new_forces

# -----------------------------
# Pygame visualization loop
# -----------------------------
pygame.init()
screen = pygame.display.set_mode((600, 600))
clock = pygame.time.Clock()

# initial forces (GPU)
forces = compute_forces(positions, velocities)

# Main loop
for step in range(steps):
    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit()
            raise SystemExit

    # integrate on GPU
    positions, velocities, forces = integrate(positions, velocities, forces)

    # render (copy to host once per frame)
    screen.fill((0, 0, 0))
    pos_np = cp.asnumpy(positions / box_size * 600.0).astype(int)
    for x, y in pos_np:
        pygame.draw.circle(screen, (255, 255, 255), (x, y), 2)
    pygame.display.flip()
    clock.tick(60)
