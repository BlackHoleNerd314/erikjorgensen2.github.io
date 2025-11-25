# n_body_simulation.py
#
# A GPU-accelerated N-body simulation script using CuPy for both Lennard-Jones
# and Newtonian gravity potentials.
#
# Lennard-Jones potential is solved using a cell-list method to achieve O(N)
# complexity for short-range interactions.
#
# Newtonian gravity is solved using a Fourier pseudospectral (Particle-Mesh)
# approach, which is efficient for long-range, periodic interactions.
#
# Requires: cupy, numpy, pygame
#
# This script is a demonstration and should be further optimized for specific
# use cases.

import cupy as cp
import numpy as np
import time
from cupyx.scipy.fft import fftn, ifftn
import pygame

# --- Simulation Parameters ---
# The number of particles. For GPU, we can handle thousands or even millions.
N = 137
# The side length of the periodic cubic box.
L = 100
# Timestep size.
dt = 1
# Total number of simulation steps.
NUM_STEPS = 1000

# Lennard-Jones parameters
epsilon_lj = 0
sigma_lj = 1
rcut_lj = 2.5 * sigma_lj  # Cutoff radius for Lennard-Jones

# Gravitational parameters
G = 0.1  # Gravitational constant
# Grid size for the Fourier transform method. Must be a power of 2.
GRID_SIZE = 256

# --- Pygame Parameters ---
# Window dimensions
SCREEN_WIDTH = 800
SCREEN_HEIGHT = 800
# Particle color
PARTICLE_COLOR = (255, 255, 255)  # White
BACKGROUND_COLOR = (0, 0, 0)  # Black
# Particle size in pixels
PARTICLE_SIZE = 1


# --- Helper Functions ---

def initialize_particles(n, box_size):
    """Initializes particles using linspace filling of a 3D grid within the simulation box.

    This places particles on a uniform lattice obtained from np.linspace. The grid size
    in each axis is ceil(n**(1/3)) so that we have at least n grid points; the extra
    points (if any) are discarded by slicing.

    Returns cupy arrays (positions, velocities, masses).
    """
    # Determine number of points per axis to fit at least n particles
    n_per_axis = int(np.ceil(n ** (1.0 / 3.0)))

    # Use endpoint=False so positions are in [0, box_size)
    coords = np.linspace(0.0, box_size, n_per_axis, endpoint=False, dtype=np.float64)
    X, Y, Z = np.meshgrid(coords, coords, coords, indexing='ij')

    positions_np = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    # Slice to first n particles (meshgrid guaranteed to have >= n points)
    positions_np = positions_np[:n].astype(np.float64)

    # Convert to cupy arrays for GPU computation
    positions = cp.asarray(positions_np)
    velocities = cp.zeros((n, 3), dtype=cp.float64)
    masses = cp.ones(n, dtype=cp.float64) / n

    return positions, velocities, masses


# A simple Velocity Verlet integrator step.
def velocity_verlet_step(positions, velocities, forces, dt):
    """
    Performs a single Velocity Verlet integration step.

    Args:
        positions (cupy.ndarray): Current particle positions.
        velocities (cupy.ndarray): Current particle velocities.
        forces (cupy.ndarray): Current forces on particles.
        dt (float): Timestep size.

    Returns:
        tuple: Updated positions and velocities.
    """
    velocities += 0.5 * forces * dt
    positions += velocities * dt
    # Apply periodic boundary conditions
    positions %= L

    return positions, velocities


# --- Lennard-Jones Force Calculation (Direct Summation Method) ---

def lennard_jones_force_direct(positions, L, rcut, epsilon, sigma):
    """
    Calculates Lennard-Jones forces using direct N^2 summation.

    Args:
        positions (cupy.ndarray): Particle positions (N, 3).
        L (float): Box side length.
        rcut (float): Cutoff radius.
        epsilon (float): LJ epsilon parameter.
        sigma (float): LJ sigma parameter.

    Returns:
        cupy.ndarray: The forces on each particle (N, 3).
    """
    forces = cp.zeros_like(positions)
    N_particles = positions.shape[0]

    # CuPy custom kernel for force calculation using direct summation
    direct_lj_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void lennard_jones_direct_kernel(
        const int N,
        const double L,
        const double rcut_sq,
        const double epsilon,
        const double sigma_6,
        const double* __restrict__ positions,
        double* __restrict__ forces) {

        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N) return;

        double3 pos_i = {positions[i*3], positions[i*3+1], positions[i*3+2]};

        double3 force_i = {0.0, 0.0, 0.0};

        // Loop over all other particles
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;

            double3 pos_j = {positions[j*3], positions[j*3+1], positions[j*3+2]};

            // Calculate separation vector with periodic boundary conditions
            double3 dr;
            dr.x = pos_i.x - pos_j.x;
            dr.y = pos_i.y - pos_j.y;
            dr.z = pos_i.z - pos_j.z;

            if (dr.x > L * 0.5) dr.x -= L;
            if (dr.x < -L * 0.5) dr.x += L;
            if (dr.y > L * 0.5) dr.y -= L;
            if (dr.y < -L * 0.5) dr.y += L;
            if (dr.z > L * 0.5) dr.z -= L;
            if (dr.z < -L * 0.5) dr.z += L;

            double r_sq = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

            if (r_sq < rcut_sq && r_sq > 1e-9) { // Avoid self-interaction and large forces at r=0
                double r_2 = 1.0 / r_sq;
                double r_6 = r_2 * r_2 * r_2;
                double r_12 = r_6 * r_6;
                double scalar_force = 24.0 * epsilon * (2.0 * sigma_6 * r_12 - sigma_6 * r_6) * r_2;

                force_i.x += scalar_force * dr.x;
                force_i.y += scalar_force * dr.y;
                force_i.z += scalar_force * dr.z;
            }
        }

        // Add the calculated force to the global force array
        // This is safe because each thread is computing its own force_i
        atomicAdd(&forces[i*3], force_i.x);
        atomicAdd(&forces[i*3+1], force_i.y);
        atomicAdd(&forces[i*3+2], force_i.z);
    }
    ''', 'lennard_jones_direct_kernel')

    threads_per_block = 256
    blocks_per_grid = (N_particles + threads_per_block - 1) // threads_per_block

    direct_lj_kernel((blocks_per_grid,), (threads_per_block,), (
        N_particles,
        L,
        rcut ** 2,
        epsilon,
        sigma ** 6,
        positions,
        forces
    ))

    return forces


# --- Lennard-Jones Force Calculation (Cell List Method) ---

def lennard_jones_force_cell_list(positions, L, rcut, epsilon, sigma):
    """
    Calculates Lennard-Jones forces using a cell list method for efficiency.

    Args:
        positions (cupy.ndarray): Particle positions (N, 3).
        L (float): Box side length.
        rcut (float): Cutoff radius.
        epsilon (float): LJ epsilon parameter.
        sigma (float): LJ sigma parameter.

    Returns:
        cupy.ndarray: The forces on each particle (N, 3).
    """
    forces = cp.zeros_like(positions)
    N_particles = positions.shape[0]

    # Define cell grid
    cell_size = rcut
    num_cells = int(np.floor(L / cell_size))

    if num_cells < 3:
        # Fallback to direct summation if cells are too large
        print("Warning: Cell size too large. Falling back to direct summation.")
        return lennard_jones_force_direct(positions, L, rcut, epsilon, sigma)

    # Map particles to cells
    cell_indices = cp.floor(positions / cell_size).astype(cp.int32)
    flat_indices = (cell_indices[:, 0] * num_cells * num_cells +
                    cell_indices[:, 1] * num_cells +
                    cell_indices[:, 2])

    # Sort particles by their cell index
    sort_perm = cp.argsort(flat_indices)
    sorted_positions = positions[sort_perm]
    sorted_cell_indices = flat_indices[sort_perm]

    # Create an index for the start and end of each cell's particle list
    cell_start_indices = cp.searchsorted(sorted_cell_indices, cp.arange(num_cells ** 3))

    # CuPy custom kernel for force calculation
    lj_kernel = cp.RawKernel(r'''
    extern "C" __global__
    void lennard_jones_cell_kernel(
        const int N,
        const double L,
        const double rcut_sq,
        const double epsilon,
        const double sigma_6,
        const int num_cells,
        const double cell_size,
        const int* __restrict__ cell_start_indices,
        const double* __restrict__ sorted_positions,
        const int* __restrict__ sort_perm,
        double* __restrict__ forces) {

        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N) return;

        // Get position of particle i
        double3 pos_i = {sorted_positions[i*3], sorted_positions[i*3+1], sorted_positions[i*3+2]};

        // Get cell index of particle i
        int cell_x = (int)floor(pos_i.x / cell_size);
        int cell_y = (int)floor(pos_i.y / cell_size);
        int cell_z = (int)floor(pos_i.z / cell_size);

        double3 force_i = {0.0, 0.0, 0.0};

        // Loop over neighboring cells (including the current one)
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {

                    int neighbor_cell_x = (cell_x + dx + num_cells) % num_cells;
                    int neighbor_cell_y = (cell_y + dy + num_cells) % num_cells;
                    int neighbor_cell_z = (cell_z + dz + num_cells) % num_cells;

                    int neighbor_cell_idx = neighbor_cell_x * num_cells * num_cells +
                                            neighbor_cell_y * num_cells +
                                            neighbor_cell_z;

                    int start_j = cell_start_indices[neighbor_cell_idx];
                    int end_j = (neighbor_cell_idx == num_cells*num_cells*num_cells - 1) ? N : cell_start_indices[neighbor_cell_idx + 1];

                    // Loop over particles in the neighbor cell
                    for (int j = start_j; j < end_j; ++j) {
                        if (i == j) continue;

                        double3 pos_j = {sorted_positions[j*3], sorted_positions[j*3+1], sorted_positions[j*3+2]};

                        // Calculate separation vector with periodic boundary conditions
                        double3 dr;
                        dr.x = pos_i.x - pos_j.x;
                        dr.y = pos_i.y - pos_j.y;
                        dr.z = pos_i.z - pos_j.z;

                        if (dr.x > L * 0.5) dr.x -= L;
                        if (dr.x < -L * 0.5) dr.x += L;
                        if (dr.y > L * 0.5) dr.y -= L;
                        if (dr.y < -L * 0.5) dr.y += L;
                        if (dr.z > L * 0.5) dr.z -= L;
                        if (dr.z < -L * 0.5) dr.z += L;

                        double r_sq = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;

                        if (r_sq < rcut_sq && r_sq > 1e-9) { // Avoid self-interaction and large forces at r=0
                            double r_2 = 1.0 / r_sq;
                            double r_6 = r_2 * r_2 * r_2;
                            double r_12 = r_6 * r_6;
                            double scalar_force = 24.0 * epsilon * (2.0 * sigma_6 * r_12 - sigma_6 * r_6) * r_2;

                            force_i.x += scalar_force * dr.x;
                            force_i.y += scalar_force * dr.y;
                            force_i.z += scalar_force * dr.z;
                        }
                    }
                }
            }
        }

        // Atomically add the calculated force to the global force array
        atomicAdd(&forces[sort_perm[i]*3], force_i.x);
        atomicAdd(&forces[sort_perm[i]*3+1], force_i.y);
        atomicAdd(&forces[sort_perm[i]*3+2], force_i.z);
    }
    ''', 'lennard_jones_cell_kernel')

    # Allocate force array and run the kernel
    forces_unsorted = cp.zeros_like(positions)
    threads_per_block = 256
    blocks_per_grid = (N_particles + threads_per_block - 1) // threads_per_block

    # Launch the kernel
    lj_kernel((blocks_per_grid,), (threads_per_block,), (
        N_particles,
        L,
        rcut ** 2,
        epsilon,
        sigma ** 6,
        num_cells,
        cell_size,
        cell_start_indices,
        sorted_positions,
        sort_perm,
        forces_unsorted
    ))

    return forces_unsorted


# --- Gravitational Force Calculation (Fourier Pseudospectral Method) ---

def gravitational_force_fourier(positions, L, G, grid_size):
    """
    Calculates gravitational forces using a Fourier pseudospectral method.

    Args:
        positions (cupy.ndarray): Particle positions (N, 3).
        L (float): Box side length.
        G (float): Gravitational constant.
        grid_size (int): Resolution of the grid.

    Returns:
        cupy.ndarray: The forces on each particle (N, 3).
    """
    # 1. Deposit particles onto a density grid
    rho_grid = cp.zeros((grid_size, grid_size, grid_size), dtype=cp.float64)
    # Simple nearest-grid-point (NGP) deposition
    indices = cp.floor(positions / L * grid_size).astype(cp.int32)
    cp.add.at(rho_grid, (indices[:, 0], indices[:, 1], indices[:, 2]), 1.0)

    # 2. Go to Fourier space
    rho_fft = fftn(rho_grid)

    # 3. Create the Green's function (gravitational kernel) in Fourier space
    k = cp.fft.fftfreq(grid_size, d=L / grid_size)
    k_x, k_y, k_z = cp.meshgrid(k, k, k, indexing='ij')
    k_sq = k_x ** 2 + k_y ** 2 + k_z ** 2

    # The kernel is 1/k^2, with special handling for k=0 (the DC component)
    kernel_fft = cp.zeros_like(k_sq)
    kernel_fft[k_sq > 0] = -4.0 * np.pi * G / (k_sq[k_sq > 0])
    kernel_fft[0, 0, 0] = 0.0  # Set DC component to zero

    # 4. Solve the Poisson equation in Fourier space
    phi_fft = rho_fft * kernel_fft

    # 5. Inverse Fourier Transform to get the potential
    phi_grid = cp.real(ifftn(phi_fft))

    # 6. Differentiate to get the force components
    grid_spacing = L / grid_size

    # Use central finite differences for the gradient
    force_grid_x = (cp.roll(phi_grid, -1, axis=0) - cp.roll(phi_grid, 1, axis=0)) / (2.0 * grid_spacing)
    force_grid_y = (cp.roll(phi_grid, -1, axis=1) - cp.roll(phi_grid, 1, axis=1)) / (2.0 * grid_spacing)
    force_grid_z = (cp.roll(phi_grid, -1, axis=2) - cp.roll(phi_grid, 1, axis=2)) / (2.0 * grid_spacing)

    # 7. Interpolate forces back to particle positions
    # We use a simple nearest-grid-point (NGP) interpolation
    forces = cp.zeros_like(positions)
    force_indices = cp.floor(positions / L * grid_size).astype(cp.int32)

    forces[:, 0] = force_grid_x[force_indices[:, 0], force_indices[:, 1], force_indices[:, 2]]
    forces[:, 1] = force_grid_y[force_indices[:, 0], force_indices[:, 1], force_indices[:, 2]]
    forces[:, 2] = force_grid_z[force_indices[:, 0], force_indices[:, 1], force_indices[:, 2]]

    # The force is the negative gradient, so we flip the sign
    return -forces


def draw_particles(positions, screen):
    """
    Draws particles onto the Pygame screen.

    Args:
        positions (cupy.ndarray): Current particle positions (N, 3).
        screen (pygame.Surface): The Pygame screen surface.
    """
    # Convert CuPy positions to a NumPy array for drawing
    positions_np = positions.get()

    screen.fill(BACKGROUND_COLOR)

    # Scale and translate positions to screen coordinates
    # We'll project the 3D particles onto the 2D xy-plane
    for i in range(positions_np.shape[0]):
        x = int(positions_np[i, 0] / L * SCREEN_WIDTH)
        y = int(positions_np[i, 1] / L * SCREEN_HEIGHT)
        pygame.draw.circle(screen, PARTICLE_COLOR, (x, y), PARTICLE_SIZE)

    pygame.display.flip()


# --- Main Simulation Loop ---

def main():
    """
    Main function to run the simulation.
    """
    print(f"Starting N-Body simulation with N={N} particles.")

    # Initialize Pygame
    pygame.init()
    screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
    pygame.display.set_caption("N-Body Simulation")
    clock = pygame.time.Clock()
    running = True

    # Initialize particle positions and velocities on the GPU using linspace-based grid
    positions, velocities, masses = initialize_particles(N, L)
    forces = cp.zeros((N, 3), dtype=cp.float64)

    # Set up which force to use for each step.
    force_type = 'LJ'  # 'LJ' or 'Gravity'

    start_time = time.time()

    for step in range(NUM_STEPS):
        print(f"Step {step + 1}/{NUM_STEPS}, Force Type: {force_type}")

        # Pygame event handling
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        if not running:
            break

        # 1. Compute forces based on the current potential
        if force_type == 'LJ':
            forces = lennard_jones_force_cell_list(positions, L, rcut_lj, epsilon_lj, sigma_lj)
            # Switch to Gravity for the next step
            force_type = 'Gravity'
        else:  # force_type == 'Gravity'
            forces = gravitational_force_fourier(positions, L, G, GRID_SIZE)
            # Switch to LJ for the next step
            force_type = 'LJ'

        # 2. Perform the Velocity Verlet integration
        positions, velocities = velocity_verlet_step(positions, velocities, forces, dt)

        # 3. Draw particles to the screen
        draw_particles(positions, screen)
        clock.tick(60)  # Limit to 60 FPS

    end_time = time.time()

    print(f"\nSimulation complete after {NUM_STEPS} steps.")
    print(f"Total time elapsed: {end_time - start_time:.4f} seconds.")
    print(f"Average time per step: {(end_time - start_time) / NUM_STEPS:.4f} seconds.")

    # Clean up Pygame
    pygame.quit()


# Entry point of the script
if __name__ == "__main__":
    try:
        main()
    except cp.cuda.runtime.CUDARuntimeError as e:
        print(f"CUDA Error: {e}. Ensure you have a compatible NVIDIA GPU and the CUDA Toolkit installed.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
