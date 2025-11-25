import cupy as cp
import numpy as np
import pygame
import sys
import math
import cupyx.scipy.fft as fft

# --- Simulation Parameters ---
GRID_SIZE = 256
NUM_PARTICLES = 4095  # Try increasing this to see the performance of the GPU
BOX_SIZE = 256  # The size of our 2D simulation box (float)
TIME_STEP = 0.05
SOFTENING_LENGTH = 1.0  # Used to prevent infinite forces in gravity calculation
MASS = 1.0
k_damping = 0.0
# Lennard-Jones parameters
SIGMA = 1.0
EPSILON = -0.1
LJ_CUTOFF = 3.0 * SIGMA  # Cutoff radius for Lennard-Jones interactions

# Gravity parameters
GRAVITATIONAL_CONSTANT = 0.001

# Pygame display parameters
SCREEN_WIDTH, SCREEN_HEIGHT = 800, 800
PARTICLE_RADIUS = 1

# --- CUDA Kernel for Lennard-Jones Force with Cell Lists ---
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

                double r2 = dx * dx + dy * dy + 1e-12; // avoid div by zero

                if (r2 < cutoff2) {
                    double r2_inv = 1.0 / r2;
                    double r6_inv = r2_inv * r2_inv * r2_inv;
                    double r12_inv = r6_inv * r6_inv;
                    double force_mag = 24.0 * epsilon * (2.0 * r12_inv - r6_inv) * r2_inv;

                    fx += force_mag * dx;
                    fy += force_mag * dy;
                    
                    
                    double fx_damp = 0.5 * k_damping * (q_vel_x - p_vel_x);
                    double fy_damp = 0.5 * k_damping * (q_vel_y - p_vel_y);
                    fx += fx_damp;
                    fy += fy_damp;
                    
                }
            }
        }
    }

    // Add computed forces to the global force array
    atomicAdd(&force[tid * 2], fx);
    atomicAdd(&force[tid * 2 + 1], fy);
}
''', 'cell_list_lj_kernel')


# --- Main Simulation Class ---
class NBodySimulation:
    def __init__(self):
        # Initialize CuPy arrays on the GPU (float64)
        self.pos = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)
        self.vel = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)  # zero velocity as requested
        self.acc = cp.zeros((NUM_PARTICLES, 2), dtype=cp.float64)
        self.mass = cp.full(NUM_PARTICLES, MASS, dtype=cp.float64)
        self.grid_size = GRID_SIZE

        # Cell list parameters
        self.cell_size_x = float(LJ_CUTOFF)
        self.cell_size_y = float(LJ_CUTOFF)
        self.num_cells_x = max(1, int(BOX_SIZE / self.cell_size_x))
        self.num_cells_y = max(1, int(BOX_SIZE / self.cell_size_y))
        self.num_cells = self.num_cells_x * self.num_cells_y

        # Pygame initialization
        pygame.init()
        self.screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))
        pygame.display.set_caption("2D N-Body Simulation (Lennard-Jones + Gravity)")
        self.clock = pygame.time.Clock()

        self.initialize_particles_on_grid()

    def initialize_particles_on_grid(self):
        """Initializes particles on a uniform linspace grid (no rounding) with zero initial velocities."""
        # Ensure we have at least NUM_PARTICLES grid points
        num_x = int(math.ceil(math.sqrt(NUM_PARTICLES)))
        num_y = int(math.ceil(NUM_PARTICLES / num_x))

        # Use linspace (float) over [0, BOX_SIZE) with no rounding
        x_coords = cp.linspace(0.0, BOX_SIZE, num_x, endpoint=False, dtype=cp.float64)
        y_coords = cp.linspace(0.0, BOX_SIZE, num_y, endpoint=False, dtype=cp.float64)

        x_grid, y_grid = cp.meshgrid(x_coords, y_coords)
        pos_flat = cp.stack([x_grid.ravel(), y_grid.ravel()], axis=1)

        # Place only the first NUM_PARTICLES positions (the rest are ignored)
        self.pos[:NUM_PARTICLES, :] = pos_flat[:NUM_PARTICLES, :]
        # Velocities are already zero; ensure dtype is float64
        self.vel[:] = 0.0

    def build_cell_list(self):
        """Builds a cell list and returns the flattened particle index list plus cell start/count arrays.
        The returned particle_cell_idx is simply the indices [0..N-1] for the *sorted* arrays so the CUDA
        kernel can index into the (already sorted) pos/vel arrays directly.
        """
        # Compute cell coordinates for each particle
        cell_x = cp.floor(self.pos[:, 0] / self.cell_size_x).astype(cp.int32)
        cell_y = cp.floor(self.pos[:, 1] / self.cell_size_y).astype(cp.int32)

        # Clip safety
        cell_x = cp.clip(cell_x, 0, self.num_cells_x - 1)
        cell_y = cp.clip(cell_y, 0, self.num_cells_y - 1)

        cell_idx = cell_y * self.num_cells_x + cell_x

        # Sort particles by cell index
        sorted_indices = cp.argsort(cell_idx)

        # Reorder the particle arrays so particles in the same cell are contiguous
        self.pos = self.pos[sorted_indices]
        self.acc = self.acc[sorted_indices]
        self.vel = self.vel[sorted_indices]
        self.mass = self.mass[sorted_indices]

        # Recompute cell indices on the sorted order
        cell_idx_sorted = cell_idx[sorted_indices]

        # Count particles per cell and compute starts
        cell_count = cp.bincount(cell_idx_sorted, minlength=self.num_cells).astype(cp.int32)
        cell_start = cp.concatenate((cp.zeros(1, dtype=cp.int32), cp.cumsum(cell_count[:-1]).astype(cp.int32)))

        # Since arrays are sorted, the particle index list for flattened cells is simply range(0, N)
        particle_cell_idx = cp.arange(NUM_PARTICLES, dtype=cp.int32)

        return particle_cell_idx, cell_start, cell_count

    def calculate_forces(self):
        """Calculates the total force (LJ + Gravity) on all particles."""
        # Reset accelerations (force per unit mass stored in acc)
        self.acc.fill(0.0)

        # --- 1. Lennard-Jones Force (short-range) using Cell Lists ---
        particle_cell_idx, cell_start, cell_count = self.build_cell_list()

        threads_per_block = 256
        blocks = (NUM_PARTICLES + threads_per_block - 1) // threads_per_block

        # Call the CUDA kernel with the correct particle_cell_idx (for sorted arrays)
        cell_list_lj_kernel((blocks,), (threads_per_block,), (
            self.pos, self.acc, particle_cell_idx, cell_start, cell_count,
            NUM_PARTICLES, float(BOX_SIZE), float(SIGMA ** 2), float(EPSILON), float(LJ_CUTOFF ** 2),
            int(self.num_cells_x), int(self.num_cells_y), float(self.cell_size_x), float(self.cell_size_y)
        ))

        # --- 2. Gravitational Force (long-range) using Pseudospectral method ---
        # a) Create a density grid using a 2D histogram
        hist, _, _ = cp.histogram2d(self.pos[:, 0], self.pos[:, 1], bins=(self.grid_size, self.grid_size),
                                    range=[[0.0, BOX_SIZE], [0.0, BOX_SIZE]], weights=self.mass)
        density_grid = cp.array(hist, dtype=cp.float64)

        # b) Compute the Green's function in Fourier space
        k_x = cp.fft.fftfreq(self.grid_size, d=BOX_SIZE / self.grid_size)
        k_y = cp.fft.fftfreq(self.grid_size, d=BOX_SIZE / self.grid_size)
        kx, ky = cp.meshgrid(k_x, k_y)
        k_squared = kx ** 2 + ky ** 2

        green_function_fft = cp.zeros((self.grid_size, self.grid_size), dtype=cp.complex128)
        nonzero = k_squared > 0
        green_function_fft[nonzero] = -1.0 / (k_squared[nonzero] * (2.0 * math.pi) ** 2)

        # c) Solve for the potential in Fourier space
        density_fft = fft.fft2(density_grid)
        potential_fft = GRAVITATIONAL_CONSTANT * density_fft * green_function_fft

        # d) Transform back to real space to get the potential grid
        potential_grid = fft.ifft2(potential_fft).real

        # e) Calculate force by taking the negative gradient of the potential grid
        force_x_grid = -cp.gradient(potential_grid, axis=1) * (self.grid_size / BOX_SIZE)
        force_y_grid = -cp.gradient(potential_grid, axis=0) * (self.grid_size / BOX_SIZE)

        # f) Interpolate forces back to particle positions (nearest-grid-point here)
        x_indices = cp.clip((self.pos[:, 0] / BOX_SIZE * self.grid_size).astype(cp.int32), 0, self.grid_size - 1)
        y_indices = cp.clip((self.pos[:, 1] / BOX_SIZE * self.grid_size).astype(cp.int32), 0, self.grid_size - 1)

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
                pos = (int(screen_pos[i, 0]) % SCREEN_WIDTH, int(screen_pos[i, 1]) % SCREEN_HEIGHT)
                pygame.draw.circle(self.screen, (255, 255, 255), pos, PARTICLE_RADIUS)

            pygame.display.flip()
            self.clock.tick(60)

        pygame.quit()
        sys.exit()


if __name__ == "__main__":
    sim = NBodySimulation()
    sim.run()
