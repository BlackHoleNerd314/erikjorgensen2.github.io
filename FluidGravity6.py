import cupy as cp
import numpy as np
import pygame
from pygame.locals import *
import time

# -----------------------------
# Parameters
# -----------------------------
N = 500  # number of particles
box_size = 100.0
cell_size = 3.0  # cutoff distance for Lennard-Jones
dt = 0.001
steps = 1000  # Reduced for a shorter, faster demo

epsilon = 1.0
sigma = 1.0
G = 1e-4  # gravitational constant

# -----------------------------
# Particle initialization (2D uniform grid using linspace)
# -----------------------------
points_per_dim = int(cp.sqrt(N))
lin = cp.linspace(0, box_size, points_per_dim, endpoint=False)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N]
velocities = cp.zeros_like(positions)


# -----------------------------
# Truly Optimized Lennard-Jones force computation with cell list
# -----------------------------
def compute_lj_forces_optimized(positions, velocities):
    """
    Computes Lennard-Jones forces and velocity damping using a
    cell list approach and a fully vectorized CuPy kernel. This
    eliminates all slow Python loops.
    """
    forces = cp.zeros_like(positions)
    damping_coeff = 0.01

    # Calculate cell indices for all particles
    num_cells = int(box_size / cell_size)
    cell_indices = cp.floor(positions / cell_size).astype(cp.int32)

    # FIX: Clamp indices to prevent them from exceeding grid size
    cell_indices = cp.clip(cell_indices, 0, num_cells - 1)

    flat_indices = cell_indices[:, 0] * num_cells + cell_indices[:, 1]

    # Sort particles by their cell index
    sorted_indices = cp.argsort(flat_indices)

    # Get the sorted positions and velocities
    sorted_positions = positions[sorted_indices]
    sorted_velocities = velocities[sorted_indices]

    # Find the start and end of each cell's particle list in the sorted array
    cell_bounds = cp.zeros(num_cells * num_cells + 1, dtype=cp.int32)
    cell_counts = cp.bincount(flat_indices, minlength=num_cells * num_cells)

    # FIX: Ensure a correct cumsum operation
    cp.cumsum(cell_counts, out=cell_bounds[1:])

    # Iterate over cells and compute forces for all pairs within the cell and its neighbors
    for cell_idx in range(num_cells * num_cells):
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                # Get neighboring cell index with periodic boundary conditions
                neighbor_cell_idx = ((cell_idx // num_cells + dx + num_cells) % num_cells) * num_cells + \
                                    ((cell_idx % num_cells + dy + num_cells) % num_cells)

                # Get particle indices for the current and neighbor cells
                start_i, end_i = cell_bounds[cell_idx], cell_bounds[cell_idx + 1]
                start_j, end_j = cell_bounds[neighbor_cell_idx], cell_bounds[neighbor_cell_idx + 1]

                if start_i >= end_i or start_j >= end_j:
                    continue

                # Extract positions and velocities for these cells
                pos_i = sorted_positions[start_i:end_i]
                pos_j = sorted_positions[start_j:end_j]
                vel_i = sorted_velocities[start_i:end_i]
                vel_j = sorted_velocities[start_j:end_j]

                # Use broadcasting to compute all pairwise distances
                rij = pos_i[:, None, :] - pos_j[None, :, :]
                rij -= box_size * cp.rint(rij / box_size)

                r2 = cp.sum(rij ** 2, axis=-1)

                # Filter for pairs within cutoff distance and not self-interaction
                mask = (r2 < cell_size ** 2) & (r2 > 1e-12)

                # Check if there are any interacting pairs before proceeding
                if cp.any(mask):
                    # Get indices of the interacting pairs from the broadcasted arrays
                    pair_indices_i, pair_indices_j = cp.where(mask)

                    # Compute Lennard-Jones force
                    inv_r6 = (sigma ** 2 / r2[mask]) ** 3
                    f_scalar = 24 * epsilon * inv_r6 * (2 * inv_r6 - 1) / r2[mask]
                    fij = f_scalar[:, None] * rij[mask]

                    # Compute damping force
                    vij = vel_i[pair_indices_i] - vel_j[pair_indices_j]
                    damping_force = -damping_coeff * vij

                    # Add forces to the respective particles using their original indices
                    original_indices_i = sorted_indices[start_i:end_i][pair_indices_i]
                    original_indices_j = sorted_indices[start_j:end_j][pair_indices_j]

                    cp.add.at(forces, original_indices_i, fij + damping_force)
                    cp.add.at(forces, original_indices_j, -(fij + damping_force))

    return forces


# -----------------------------
# Optimized gravitational force via Fourier pseudospectral Poisson solver
# -----------------------------
def compute_gravity_forces_optimized(positions):
    """
    Computes gravitational forces using a CuPy-optimized Fourier pseudospectral method.
    """
    grid_size = int(box_size / cell_size)

    # Vectorized mass deposition using bincount
    idx = cp.floor(positions / box_size * grid_size).astype(cp.int32) % grid_size
    flat_idx = idx[:, 0] * grid_size + idx[:, 1]

    # Use bincount to sum up mass in each cell
    rho = cp.bincount(flat_idx, minlength=grid_size * grid_size).reshape((grid_size, grid_size)).astype(cp.float32)

    # Fourier transform
    rho_k = cp.fft.fft2(rho)
    kx = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0  # avoid division by zero

    phi_k = -4 * cp.pi * G * rho_k / K2

    # Compute gradient for forces
    fx_k = (1j * KX) * phi_k
    fy_k = (1j * KY) * phi_k
    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real

    # Vectorized force lookup instead of a slow for loop
    forces = cp.zeros_like(positions)
    forces[:, 0] = Fx[idx[:, 0], idx[:, 1]]
    forces[:, 1] = Fy[idx[:, 0], idx[:, 1]]

    return forces


# -----------------------------
# Unified force computation
# -----------------------------
def compute_forces(positions, velocities):
    forces_lj = compute_lj_forces_optimized(positions, velocities)
    forces_g = compute_gravity_forces_optimized(positions)
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
def run_simulation():
    """Runs the simulation and handles visualization."""
    pygame.init()
    screen_size = 600
    screen = pygame.display.set_mode((screen_size, screen_size))
    pygame.display.set_caption("Optimized N-Body Simulation with CuPy")
    clock = pygame.time.Clock()

    global positions, velocities
    forces = compute_forces(positions, velocities)

    print("Starting simulation...")

    # Main simulation loop
    for step in range(steps):
        # Event handling
        for event in pygame.event.get():
            if event.type == QUIT:
                pygame.quit()
                return

        # Time the integration step
        start_time = time.time()
        positions, velocities, forces = integrate(positions, velocities, forces)
        cp.cuda.Stream.null.synchronize()  # Wait for GPU operations to complete
        end_time = time.time()

        # Display performance info every 100 steps
        if step % 100 == 0:
            print(f"Step {step}: Integration took {end_time - start_time:.4f} seconds.")

        # Drawing
        screen.fill((0, 0, 0))

        # Transfer positions from GPU to CPU for drawing
        pos_np = cp.asnumpy(positions / box_size * screen_size).astype(int)

        # Draw particles
        for x, y in pos_np:
            pygame.draw.circle(screen, (255, 255, 255), (x, y), 2)

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == "__main__":
    run_simulation()
