# FluidGravity_safe_damped.py
# GPU fluid+gravity simulation with safety guards and tunable damping

import cupy as cp
import numpy as np
import pygame
import sys
import traceback
from pygame.locals import *

# -----------------------------
# Parameters
# -----------------------------
N_target = 4096
box_size = 1024
cell_size = 4.0
num_cells = int(box_size / cell_size)
dt = 0.005
steps = 200000

epsilon = 1.0
sigma = 1.0
G = 0.001
damping = 0.01

# Safety tuning
v_max = 10.0
eps_safe = 1e-12
save_interval = 1e10
max_ppc = 64
threads_per_block = 128

# -----------------------------
# Initialization
# -----------------------------
points_per_dim = int(cp.sqrt(N_target))
N = points_per_dim * points_per_dim
print(f"Initializing N={N} particles ({points_per_dim}x{points_per_dim} grid).")

lin = cp.linspace(0, box_size, points_per_dim, endpoint=False, dtype=cp.float32)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N].astype(cp.float32)
positions += (cp.random.rand(N, 2).astype(cp.float32) - 0.5) * 0.1
velocities = (cp.random.rand(N, 2).astype(cp.float32) - 0.5) * 10

# -----------------------------
# CUDA Kernels
# -----------------------------
neighbor_build_src = r'''
extern "C" __global__
void build_neighbor_list(const float* pos, const int N,
                         const float box_size, const float cell_size,
                         const int num_cells,
                         int* cell_particles, int* cell_counts,
                         const int max_ppc) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N) return;

    float px = pos[2*i];
    float py = pos[2*i + 1];

    int cx = (int)floorf(px / cell_size);
    int cy = (int)floorf(py / cell_size);
    cx = cx % num_cells; if (cx < 0) cx += num_cells;
    cy = cy % num_cells; if (cy < 0) cy += num_cells;

    int cidx = cx * num_cells + cy;
    int offset = atomicAdd(&cell_counts[cidx], 1);
    if (offset < max_ppc)
        cell_particles[cidx * max_ppc + offset] = i;
}
'''

force_compute_src = r'''
extern "C" __global__
void compute_local_forces(const float* pos, const float* vel, float* forces,
                          const int N, const float box_size, const float cell_size,
                          const int num_cells,
                          const int* cell_particles, const int* cell_counts,
                          const int max_ppc, const float epsilon, const float sigma,
                          const float damping, const float cutoff2, const float eps_safe) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N) return;

    float xi = pos[2*i];
    float yi = pos[2*i + 1];
    float vix = vel[2*i];
    float viy = vel[2*i + 1];

    int cx = (int)floorf(xi / cell_size);
    int cy = (int)floorf(yi / cell_size);
    cx = cx % num_cells; if (cx < 0) cx += num_cells;
    cy = cy % num_cells; if (cy < 0) cy += num_cells;

    float fx = 0.0f;
    float fy = 0.0f;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            int ncx = cx + dx;
            int ncy = cy + dy;
            if (ncx >= num_cells) ncx -= num_cells;
            if (ncx < 0) ncx += num_cells;
            if (ncy >= num_cells) ncy -= num_cells;
            if (ncy < 0) ncy += num_cells;
            int cidx = ncx * num_cells + ncy;
            int count = cell_counts[cidx];
            if (count <= 0) continue;

            int base = cidx * max_ppc;
            for (int k = 0; k < count; ++k) {
                int j = cell_particles[base + k];
                if (j < 0 || j == i) continue;

                float xj = pos[2*j];
                float yj = pos[2*j + 1];
                float vjx = vel[2*j];
                float vjy = vel[2*j + 1];

                float dxij = xi - xj;
                float dyij = yi - yj;

                float half = 0.5f * box_size;
                if (dxij > half) dxij -= box_size;
                if (dxij < -half) dxij += box_size;
                if (dyij > half) dyij -= box_size;
                if (dyij < -half) dyij += box_size;

                float r2 = dxij*dxij + dyij*dyij;
                if (r2 <= eps_safe) continue;

                if (r2 < cutoff2) {
                    float inv_r2 = 1.0f / (r2 + eps_safe);
                    float inv_r6 = (sigma * sigma * inv_r2);
                    inv_r6 = inv_r6 * inv_r6 * inv_r6;
                    float f_scalar = 24.0f * epsilon * inv_r6 * (2.0f * inv_r6 - 1.0f) * inv_r2;

                    float fijx = f_scalar * dxij;
                    float fijy = f_scalar * dyij;

                    float dvx = vix - vjx;
                    float dvy = viy - vjy;
                    float dfx = -damping * dvx;
                    float dfy = -damping * dvy;

                    float total_x = fijx + dfx;
                    float total_y = fijy + dfy;

                    atomicAdd(&forces[2*i], total_x);
                    atomicAdd(&forces[2*i + 1], total_y);
                    atomicAdd(&forces[2*j], -total_x);
                    atomicAdd(&forces[2*j + 1], -total_y);
                }
            }
        }
    }
}
'''

build_neighbor_list_kernel = cp.RawKernel(neighbor_build_src, 'build_neighbor_list')
compute_local_forces_kernel = cp.RawKernel(force_compute_src, 'compute_local_forces')

# -----------------------------
# Gravity solver
# -----------------------------
def compute_gravity_forces(positions):
    grid_size = num_cells
    rho = cp.zeros((grid_size, grid_size), dtype=cp.float32)
    safe_pos = cp.nan_to_num(positions, nan=0.0, posinf=0.0, neginf=0.0)
    idx = cp.floor(safe_pos / box_size * grid_size).astype(cp.int32) % grid_size
    cp.add.at(rho, (idx[:, 0], idx[:, 1]), 1.0)

    rho_k = cp.fft.fft2(rho)
    dk = box_size / grid_size
    kx = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0
    phi_k = -4.0 * cp.pi * G * rho_k / (K2 + 1e-20)
    phi_k[0, 0] = 0.0
    fx_k = (1j * KX) * phi_k
    fy_k = (1j * KY) * phi_k
    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real
    forces = cp.zeros_like(positions, dtype=cp.float32)
    gx, gy = idx[:, 0], idx[:, 1]
    forces[:, 0] = Fx[gx, gy]
    forces[:, 1] = Fy[gx, gy]
    return -forces

# -----------------------------
# Force computation wrapper
# -----------------------------
def compute_forces(positions, velocities):
    Nloc = positions.shape[0]
    ncell = num_cells * num_cells
    cell_particles = cp.full((ncell * max_ppc,), -1, dtype=cp.int32)
    cell_counts = cp.zeros((ncell,), dtype=cp.int32)
    blocks = (Nloc + threads_per_block - 1) // threads_per_block

    build_neighbor_list_kernel((blocks,), (threads_per_block,),
        (positions.ravel(), cp.int32(Nloc), cp.float32(box_size),
         cp.float32(cell_size), cp.int32(num_cells),
         cell_particles, cell_counts, cp.int32(max_ppc)))

    forces = cp.zeros((Nloc, 2), dtype=cp.float32)
    forces_flat = forces.ravel()
    cutoff2 = float(cell_size * cell_size)

    compute_local_forces_kernel((blocks,), (threads_per_block,),
        (positions.ravel(), velocities.ravel(), forces_flat,
         cp.int32(Nloc), cp.float32(box_size), cp.float32(cell_size),
         cp.int32(num_cells), cell_particles, cell_counts,
         cp.int32(max_ppc), cp.float32(epsilon), cp.float32(sigma),
         cp.float32(damping), cp.float32(cutoff2), cp.float32(eps_safe)))

    forces += compute_gravity_forces(positions)
    return forces

# -----------------------------
# Integration & safety
# -----------------------------
def integrate(positions, velocities, forces):
    positions = (positions + velocities * dt + 0.5 * forces * dt * dt) % box_size
    positions = cp.nan_to_num(positions, nan=0.0, posinf=0.0, neginf=0.0)
    new_forces = compute_forces(positions, velocities)
    velocities = velocities + 0.5 * (forces + new_forces) * dt
    velocities = cp.nan_to_num(velocities, nan=0.0, posinf=v_max, neginf=-v_max)
    vnorm = cp.linalg.norm(velocities, axis=1, keepdims=True)
    velocities = cp.where(vnorm > v_max, velocities * (v_max / (vnorm + 1e-20)), velocities)
    return positions, velocities, new_forces

# -----------------------------
# Save snapshot
# -----------------------------
def save_snapshot(positions, velocities, step, fname):
    cp.savez(fname, positions=positions, velocities=velocities, step=np.int64(step))
    print(f"Saved {fname}")

# -----------------------------
# Main loop with Pygame visualization
# -----------------------------
def main_loop():
    pygame.init()
    screen = pygame.display.set_mode((600, 600))
    clock = pygame.time.Clock()
    forces = compute_forces(positions, velocities)

    for step in range(steps):
        for event in pygame.event.get():
            if event.type == QUIT:
                save_snapshot(positions, velocities, step, "exit_snapshot.npz")
                pygame.quit()
                return

        positions_local, velocities_local, forces_local = integrate(positions, velocities, forces)
        positions[:] = positions_local
        velocities[:] = velocities_local
        forces[:] = forces_local

        if step % 10 == 0:
            screen.fill((0, 0, 0))
            pos_np = cp.asnumpy(positions / box_size * 600).astype(int)
            pos_np = np.clip(pos_np, 0, 599)

            # Grayscale histogram
            hist = np.zeros((600, 600), dtype=np.float32)
            np.add.at(hist, (pos_np[:, 1], pos_np[:, 0]), 1)  # Y first, then X
            hist = np.clip(hist / hist.max() * 255.0, 0, 255).astype(np.uint8)

            # Convert to surface
            surf = pygame.surfarray.make_surface(np.stack([hist] * 3, axis=-1))
            screen.blit(surf, (0, 0))
            pygame.display.flip()
            clock.tick(12)

        if step % save_interval == 0 and step > 0:
            save_snapshot(positions, velocities, step, f"snapshot_{step}.npz")

    save_snapshot(positions, velocities, steps, "final_snapshot.npz")
    pygame.quit()

if __name__ == "__main__":
    main_loop()
