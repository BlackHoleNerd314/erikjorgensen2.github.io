import cupy as cp
import numpy as np
import pygame
from pygame.locals import *

# -----------------------------
# Parameters
# -----------------------------
N = 1024  # number of particles (32x32 grid)
box_size = 100.0
cell_size = 3.0  # cutoff distance for Lennard-Jones
num_cells = int(box_size / cell_size)
dt = 0.001
steps = 100000

epsilon = 1.0
sigma = 1.0
G = 0.0001  # gravitational constant

# -----------------------------
# Particle initialization (uniform grid + random velocities)
# -----------------------------
points_per_dim = int(cp.sqrt(N))
lin = cp.linspace(0, box_size, points_per_dim, endpoint=False)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N]
velocities = (cp.random.rand(N, 2) - 0.5) * 2.0  # random initial velocity

# -----------------------------
# GPU kernel: build neighbor list (spatial hashing)
# -----------------------------
neigh_kernel = cp.RawKernel(r'''
extern "C" __global__ void build_neighbor_list(
    const float* pos, const int N, const float box_size,
    const float cell_size, const int num_cells,
    int* cell_particles, int* cell_counts, int max_ppc) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N) return;
    int cx = (int)floorf(pos[2*i] / cell_size);
    int cy = (int)floorf(pos[2*i+1] / cell_size);
    cx = (cx % num_cells + num_cells) % num_cells;
    cy = (cy % num_cells + num_cells) % num_cells;
    int cell_index = cx * num_cells + cy;
    int offset = atomicAdd(&cell_counts[cell_index], 1);
    if (offset < max_ppc) {
        cell_particles[cell_index * max_ppc + offset] = i;
    }
}
''', 'build_neighbor_list')

# -----------------------------
# Lennard-Jones + short-range damping kernel
# -----------------------------
force_kernel = cp.RawKernel(r'''
extern "C" __global__ void compute_forces(
    const float* pos, const float* vel, float* forces,
    const int N, const float box_size, const float cell_size,
    const int num_cells, const int* cell_particles, const int* cell_counts,
    const int max_ppc, const float epsilon, const float sigma, const float damping) {

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N) return;
    float fx = 0.0f, fy = 0.0f;

    int cx = (int)floorf(pos[2*i] / cell_size);
    int cy = (int)floorf(pos[2*i+1] / cell_size);

    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            int ncx = (cx + dx + num_cells) % num_cells;
            int ncy = (cy + dy + num_cells) % num_cells;
            int cell_index = ncx * num_cells + ncy;
            int count = cell_counts[cell_index];

            for (int k = 0; k < count; k++) {
                int j = cell_particles[cell_index * max_ppc + k];
                if (j == i) continue;
                float rx = pos[2*i] - pos[2*j];
                float ry = pos[2*i+1] - pos[2*j+1];
                rx -= box_size * roundf(rx / box_size);
                ry -= box_size * roundf(ry / box_size);

                float r2 = rx*rx + ry*ry;
                if (r2 > 1e-12f && r2 < cell_size*cell_size) {
                    float inv_r2 = 1.0f / r2;
                    float inv_r6 = powf(sigma*sigma * inv_r2, 3);
                    float f_scalar = 24.0f * epsilon * inv_r6 * (2*inv_r6 - 1) * inv_r2;
                    fx += f_scalar * rx;
                    fy += f_scalar * ry;

                    float vx = vel[2*i] - vel[2*j];
                    float vy = vel[2*i+1] - vel[2*j+1];
                    fx -= damping * vx;
                    fy -= damping * vy;
                }
            }
        }
    }
    forces[2*i] = fx;
    forces[2*i+1] = fy;
}
''', 'compute_forces')


# -----------------------------
# Gravitational force via FFT
# -----------------------------
def compute_gravity_forces(positions):
    grid_size = num_cells
    rho = cp.zeros((grid_size, grid_size), dtype=cp.float32)
    idx = cp.floor(positions / box_size * grid_size).astype(cp.int32) % grid_size
    for i in range(len(positions)):
        rho[idx[i, 0], idx[i, 1]] += 1.0

    rho_k = cp.fft.fft2(rho)
    kx = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=box_size / grid_size) * 2 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')
    K2 = KX ** 2 + KY ** 2
    K2[0, 0] = 1.0

    phi_k = -4 * cp.pi * G * rho_k / K2
    fx_k = 1j * KX * phi_k
    fy_k = 1j * KY * phi_k

    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real

    forces = cp.zeros_like(positions)
    for i in range(len(positions)):
        gx, gy = idx[i]
        forces[i, 0] = Fx[gx, gy]
        forces[i, 1] = Fy[gx, gy]
    return forces


# -----------------------------
# Unified force computation (GPU)
# -----------------------------
def compute_forces(positions, velocities):
    N = len(positions)
    threads = 128
    blocks = (N + threads - 1) // threads

    max_ppc = 32
    cell_particles = cp.full((num_cells * num_cells * max_ppc,), -1, dtype=cp.int32)
    cell_counts = cp.zeros((num_cells * num_cells,), dtype=cp.int32)

    neigh_kernel((blocks,), (threads,),
                 (positions.astype(cp.float32).ravel(), N, box_size, cell_size, num_cells,
                  cell_particles, cell_counts, max_ppc))

    forces = cp.zeros_like(positions, dtype=cp.float32)
    force_kernel((blocks,), (threads,),
                 (positions.astype(cp.float32).ravel(), velocities.astype(cp.float32).ravel(),
                  forces.ravel(), N, box_size, cell_size, num_cells, cell_particles, cell_counts,
                  max_ppc, epsilon, sigma, 0.01))

    forces += compute_gravity_forces(positions)
    return forces


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
# Visualization (pygame)
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

    if step % 5 == 0:
        screen.fill((0, 0, 0))
        pos_np = cp.asnumpy(positions / box_size * 600).astype(int)
        for x, y in pos_np:
            pygame.draw.circle(screen, (255, 255, 255), (x, y), 2)
        pygame.display.flip()
        clock.tick(60)

# -----------------------------
# Save snapshot to file
# -----------------------------
cp.savez('fluid_gravity_gpu_snapshot.npz', positions=positions, velocities=velocities)
