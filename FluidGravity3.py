import cupy as cp
import numpy as np
import pygame
import sys
import math
import cupyx.scipy.fft as fft

# --- Simulation Parameters ---
# Use a power of 2 for FFT grid size for efficiency.
GRID_SIZE = 256
NUM_PARTICLES = 512  # Try increasing this to see the performance of the GPU
BOX_SIZE = 4096  # The size of our 2D simulation box
TIME_STEP = 0.005
SOFTENING_LENGTH = 1  # Used to prevent infinite forces in gravity calculation
MASS = 1.0

# Lennard-Jones parameters
SIGMA = 1
EPSILON = 1
LJ_CUTOFF = 3.0 * SIGMA  # Cutoff radius for Lennard-Jones interactions

# Gravity parameters
GRAVITATIONAL_CONSTANT = 5

# Pygame display parameters
SCREEN_WIDTH, SCREEN_HEIGHT = 800, 800
PARTICLE_RADIUS = 2

# --- CUDA Kernel for Lennard-Jones Force with Cell Lists ---
# This kernel calculates Lennard-Jones force by iterating over a particle's neighbors
# identified by a cell list. This is much more efficient than the brute-force method.
cell_list_lj_kernel = cp.RawKernel(r'''
extern "C" __global__
void cell_list_lj_kernel(
    const double* pos, 
    double* force, 
    const int* particle_cell_idx,
    const int* cell_start,
    const int* cell_count,
    const int N, 
    const double box_size, 
    const double sigma2, 
    const double epsilon,
    const double cutoff2,
    const int num_cells_x,
    const int num_cells_y,
    const double cell_size_x,
    const double cell_size_y) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= N) return;

    double p_x = pos[tid * 2];
    double p_y = pos[tid * 2 + 1];

    int p_cell_x = (int)(p_x / cell_size_x);
    int p_cell_y = (int)(p_y / cell_size_y);

    double fx = 0.0;
    double fy = 0.0;

    // Iterate over neighboring cells (including the particle's own cell)
    for (int i_x = -1; i_x <= 1; ++i_x) {
        for (int i_y = -1; i_y <= 1; ++i_y) {
            int neighbor_cell_x = (p_cell_x + i_x + num_cells_x) % num_cells_x;
            int neighbor_cell_y = (p_cell_y + i_y + num_cells_y) % num_cells_y;

            int neighbor_cell_flat_idx = neighbor_cell_y * num_cells_x + neighbor_cell_x;

            int start_idx = cell_start[neighbor_cell_flat_idx];
            int end_idx = start_idx + cell_count[neighbor_cell_flat_idx];

            for (int j = start_idx; j < end_idx; ++j) {
                int neighbor_tid = particle_cell_idx[j];
                if (tid == neighbor_tid) continue;

                double q_x = pos[neighbor_tid * 2];
                double q_y = pos[neighbor_tid * 2 + 1];

                // Periodic boundary conditions
                double dx = q_x - p_x;
                double dy = q_y - p_y;

                if (dx > box_size / 2.0) dx -= box_size;
                if (dx < -box_size / 2.0) dx += box_size;
                if (dy > box_size / 2.0) dy -= box_size;
                if (dy < -box_size / 2.0) dy += box_size;

                double r2 = dx * dx + dy * dy;

                if (r2 < cutoff2) {
                    double r2_inv = 1.0 / r2;
                    double r6_inv = r2_inv * r2_inv * r2_inv;
                    double r12_inv = r6_inv * r6_inv;
                    double force_mag = 24.0 * epsilon * (2.0 * r12_inv - r6_inv) * r2_inv;

                    fx += force_mag * dx;
                    fy += force_mag * dy;
                }
            }
        }
    }

    // Add computed forces to the global force array
    force[tid * 2] += fx;
    force[tid * 2 + 1] += fy;
}
''', 'cell_list_lj_kernel')


# --- Main Simulation Class ---
class NBodySimulation:
    def __init__(self):
        # Initialize CuPy arrays on the GPU
        self.pos = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)
        self.vel = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)
        self.acc = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)
        self.mass = cp.full(NUM_PARTICLES, MASS, dtype=cp.float64)
        self.grid_size = GRID_SIZE

        # Cell list parameters
        self.cell_size_x = LJ_CUTOFF
        self.cell_size_y = LJ_CUTOFF
        self.num_cells_x = int(BOX_SIZE / self.cell_size_x)
        self.num_cells_y = int(BOX_SIZE / self.cell_size_y)
        self.num_cells = self.num_cells_x * self.num_cells_y

        # Pygame initialization
        pygame.init()
        self.screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
        pygame.display.set_caption("2D N-Body Simulation (Lennard-Jones + Gravity)")
        self.clock = pygame.time.Clock()

        self.initialize_particles_on_grid()

    def initialize_particles_on_grid(self):
        """Initializes particles on a uniform grid to avoid initial close-range forces."""
        num_x = int(math.sqrt(NUM_PARTICLES))
        num_y = int(NUM_PARTICLES / num_x)

        x_coords = cp.linspace(0, BOX_SIZE, num_x, endpoint=False)
        y_coords = cp.linspace(0, BOX_SIZE, num_y, endpoint=False)

        x_grid, y_grid = cp.meshgrid(x_coords, y_coords)

        # Reshape and set positions
        pos_flat = cp.stack([x_grid.ravel(), y_grid.ravel()], axis=1)
        self.pos[:pos_flat.shape[0], :] = pos_flat

    def build_cell_list(self):
        """Builds a cell list to group particles by spatial location for efficient neighbor finding."""

        # Get cell indices for each particle
        cell_x = (self.pos[:, 0] / self.cell_size_x).astype(cp.int32)
        cell_y = (self.pos[:, 1] / self.cell_size_y).astype(cp.int32)

        # Flatten the 2D cell index to a 1D index
        cell_idx = cell_y * self.num_cells_x + cell_x

        # Sort particles by their cell index
        sorted_indices = cp.argsort(cell_idx)
        self.pos = self.pos[sorted_indices]
        self.acc = self.acc[sorted_indices]
        self.vel = self.vel[sorted_indices]
        self.mass = self.mass[sorted_indices]

        # Re-calculate the cell indices for the sorted particles
        cell_idx_sorted = cell_idx[sorted_indices]

        # Compute the number of particles in each cell
        cell_count = cp.bincount(cell_idx_sorted, minlength=self.num_cells)

        # Compute the starting index for each cell
        cell_start = cp.concatenate((cp.zeros(1, dtype=cp.int32), cp.cumsum(cell_count)[:-1]))

        return cell_idx_sorted, cell_start, cell_count

    def calculate_forces(self):
        """Calculates the total force (LJ + Gravity) on all particles."""
        self.acc.fill(0.0)  # Reset accelerations

        # --- 1. Lennard-Jones Force (short-range) using Cell Lists ---
        # Build the cell list on the GPU
        particle_cell_idx, cell_start, cell_count = self.build_cell_list()

        # Call the new CUDA kernel with cell list data
        threads_per_block = 256
        blocks = (NUM_PARTICLES + threads_per_block - 1) // threads_per_block
        cell_list_lj_kernel((blocks,), (threads_per_block,), (
        self.pos, self.acc, cp.arange(NUM_PARTICLES), cell_start, cell_count, NUM_PARTICLES, BOX_SIZE, SIGMA ** 2,
        EPSILON, LJ_CUTOFF ** 2, self.num_cells_x, self.num_cells_y, self.cell_size_x, self.cell_size_y))

        # --- 2. Gravitational Force (long-range) using Pseudospectral method ---
        # Note: The positions are now sorted. We need to handle this for visualization later

        # a) Create a density grid using a 2D histogram
        hist, _, _ = cp.histogram2d(self.pos[:, 0], self.pos[:, 1], bins=(self.grid_size, self.grid_size),
                                    range=[[0, BOX_SIZE], [0, BOX_SIZE]], weights=self.mass)
        density_grid = cp.array(hist)

        # b) Compute the Green's function in Fourier space
        k_x = cp.fft.fftfreq(self.grid_size, d=BOX_SIZE / self.grid_size)
        k_y = cp.fft.fftfreq(self.grid_size, d=BOX_SIZE / self.grid_size)
        kx, ky = cp.meshgrid(k_x, k_y)
        k_squared = kx ** 2 + ky ** 2

        # Handle the zero frequency component to avoid division by zero
        green_function_fft = cp.zeros((self.grid_size, self.grid_size), dtype=cp.complex128)
        nonzero_indices = k_squared > 0
        green_function_fft[nonzero_indices] = -1.0 / (k_squared[nonzero_indices] * 4 * math.pi ** 2)

        # c) Solve for the potential in Fourier space
        density_fft = fft.fft2(density_grid)
        potential_fft = GRAVITATIONAL_CONSTANT * density_fft * green_function_fft

        # d) Transform back to real space to get the potential grid
        potential_grid = fft.ifft2(potential_fft).real

        # e) Calculate force by taking the negative gradient of the potential grid
        force_x_grid = -cp.gradient(potential_grid, axis=1) * self.grid_size / BOX_SIZE
        force_y_grid = -cp.gradient(potential_grid, axis=0) * self.grid_size / BOX_SIZE

        # f) Interpolate forces back to particle positions
        x_indices = cp.clip((self.pos[:, 0] / BOX_SIZE * self.grid_size).astype(cp.int32), 0, self.grid_size - 1)
        y_indices = cp.clip((self.pos[:, 1] / BOX_SIZE * self.grid_size).astype(cp.int32), 0, self.grid_size - 1)

        # Update the acceleration array with gravitational force
        self.acc[:, 0] += force_x_grid[y_indices, x_indices]
        self.acc[:, 1] += force_y_grid[y_indices, x_indices]

    def update_particles(self):
        """Updates particle positions and velocities using a simple Euler integration."""
        self.vel += self.acc * TIME_STEP
        self.pos += self.vel * TIME_STEP

        # Apply periodic boundary conditions
        self.pos %= BOX_SIZE

    def run(self):
        """Main simulation loop."""
        running = True
        while running:
            # --- Event Handling ---
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False

            # --- Simulation Step (on GPU) ---
            self.calculate_forces()
            self.update_particles()

            # --- Visualization (on CPU) ---
            self.screen.fill((0, 0, 0))  # Black background

            # Get particle positions from GPU to CPU
            cpu_pos = self.pos.get()

            # Map simulation coordinates to screen coordinates
            screen_pos = (cpu_pos / BOX_SIZE) * np.array([SCREEN_WIDTH, SCREEN_HEIGHT])

            # Draw particles
            for i in range(NUM_PARTICLES):
                pos = (int(screen_pos[i, 0]), int(screen_pos[i, 1]))
                pygame.draw.circle(self.screen, (255, 255, 255), pos, PARTICLE_RADIUS)

            pygame.display.flip()
            self.clock.tick(60)

            # Print FPS for performance monitoring
            # print(f"FPS: {self.clock.get_fps():.2f}")

        pygame.quit()
        sys.exit()


if __name__ == "__main__":
    sim = NBodySimulation()
    sim.run()
