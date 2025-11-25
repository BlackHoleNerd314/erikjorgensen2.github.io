# FluidGravity_safe.py
# GPU fluid+gravity simulation with safety guards to avoid illegal memory access.

import cupy as cp
import numpy as np
import pygame
import sys
import traceback
from pygame.locals import *

# -----------------------------
# Parameters
# -----------------------------
# NOTE: choose N as a perfect square for a uniform grid initialization, or N will be reduced to
# points_per_dim**2 automatically.
N_target = 4096  # requested number of particles (will become points_per_dim**2)
box_size = 100.0
cell_size = 3.0      # LJ cutoff and cell sizing
num_cells = int(box_size / cell_size)
dt = 0.001
steps = 200000

epsilon = 1.0
sigma = 1.0
G = 0.001

# Safety tuning
v_max = 100.0            # max particle speed allowed (caps runaway)
eps_safe = 1e-12         # used to protect divisions
save_interval = 100      # save snapshot every this many steps
max_ppc = 64             # max particles per cell (increase to avoid overflow)
threads_per_block = 128

# -----------------------------
# Init particles (grid + random velocity)
# -----------------------------
points_per_dim = int(cp.sqrt(N_target))
N = points_per_dim * points_per_dim
print(f"Initializing N={N} particles ({points_per_dim}x{points_per_dim} grid).")

lin = cp.linspace(0, box_size, points_per_dim, endpoint=False, dtype=cp.float32)
X, Y = cp.meshgrid(lin, lin, indexing='ij')
positions = cp.vstack([X.ravel(), Y.ravel()]).T[:N].astype(cp.float32)  # (N,2)
# small random perturbation to positions to avoid exact grid degeneracy
positions += (cp.random.rand(N, 2).astype(cp.float32) - 0.5) * 1e-3

# random initial velocities (GPU)
velocities = (cp.random.rand(N, 2).astype(cp.float32) - 0.5) * 2.0

# -----------------------------
# Raw kernels (neighbor cell deposit & local LJ/damping forces)
# The kernels assume:
#  - pos: float32[N*2] flattened (x0,y0,x1,y1,...)
#  - vel: float32[N*2]
#  - cell_particles: int32[num_cells*num_cells*max_ppc] (initialized -1)
#  - cell_counts: int32[num_cells*num_cells] (initialized 0)
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

    // load particle position
    float px = pos[2*i];
    float py = pos[2*i + 1];

    int cx = (int)floorf(px / cell_size);
    int cy = (int)floorf(py / cell_size);

    // wrap indices to [0, num_cells-1]
    cx = cx % num_cells; if (cx < 0) cx += num_cells;
    cy = cy % num_cells; if (cy < 0) cy += num_cells;

    int cidx = cx * num_cells + cy;
    int offset = atomicAdd(&cell_counts[cidx], 1);

    if (offset < max_ppc) {
        cell_particles[cidx * max_ppc + offset] = i;
    } else {
        // overflow: we silently drop additional particles in that cell.
        // This is safer than writing out-of-bounds.
    }
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

    // iterate 3x3 neighbor cells
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

                // minimum image
                float half = 0.5f * box_size;
                if (dxij > half) dxij -= box_size;
                if (dxij < -half) dxij += box_size;
                if (dyij > half) dyij -= box_size;
                if (dyij < -half) dyij += box_size;

                float r2 = dxij*dxij + dyij*dyij;
                if (r2 <= eps_safe) continue; // skip too-close or zero

                if (r2 < cutoff2) {
                    float inv_r2 = 1.0f / (r2 + eps_safe);
                    float inv_r6 = (sigma * sigma * inv_r2);
                    inv_r6 = inv_r6 * inv_r6 * inv_r6; // (sigma^6/r^6)
                    float f_scalar = 24.0f * epsilon * inv_r6 * (2.0f * inv_r6 - 1.0f) * inv_r2;

                    float fijx = f_scalar * dxij;
                    float fijy = f_scalar * dyij;

                    float dvx = vix - vjx;
                    float dvy = viy - vjy;
                    float dfx = -damping * dvx;
                    float dfy = -damping * dvy;

                    float total_x = fijx + dfx;
                    float total_y = fijy + dfy;

                    // Accumulate to forces[i] atomically to avoid race with other threads updating same i.
                    atomicAdd(&forces[2*i], total_x);
                    atomicAdd(&forces[2*i + 1], total_y);

                    // Subtract from j atomically (both threads may handle pair, but safe)
                    atomicAdd(&forces[2*j], -total_x);
                    atomicAdd(&forces[2*j + 1], -total_y);
                }
            } // for k
        }
    } // neighbor loops
}
'''

# compile kernels
build_neighbor_list_kernel = cp.RawKernel(neighbor_build_src, 'build_neighbor_list')
compute_local_forces_kernel = cp.RawKernel(force_compute_src, 'compute_local_forces')

# -----------------------------
# Gravity solver (FFT-based). Uses safe indexing and nan-handling.
# -----------------------------
def compute_gravity_forces(positions):
    grid_size = num_cells
    rho = cp.zeros((grid_size, grid_size), dtype=cp.float32)

    # safe positions -> indices (guard NaN/inf)
    safe_pos = cp.nan_to_num(positions, nan=0.0, posinf=0.0, neginf=0.0)
    idx = cp.floor(safe_pos / box_size * grid_size).astype(cp.int32) % grid_size

    # deposit mass using cp.add.at (GPU)
    cp.add.at(rho, (idx[:, 0], idx[:, 1]), 1.0)

    rho_k = cp.fft.fft2(rho)
    dk = box_size / grid_size
    kx = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    ky = cp.fft.fftfreq(grid_size, d=dk) * 2.0 * cp.pi
    KX, KY = cp.meshgrid(kx, ky, indexing='ij')
    K2 = KX**2 + KY**2
    K2[0, 0] = 1.0  # avoid division by zero, we'll zero phi_k[0,0] below

    phi_k = -4.0 * cp.pi * G * rho_k / (K2 + 1e-20)
    phi_k[0, 0] = 0.0

    fx_k = (1j * KX) * phi_k
    fy_k = (1j * KY) * phi_k

    Fx = cp.fft.ifft2(fx_k).real
    Fy = cp.fft.ifft2(fy_k).real

    forces = cp.zeros_like(positions, dtype=cp.float32)
    # Use safe_pos for indexing to avoid NaNs
    gx = idx[:, 0] % grid_size
    gy = idx[:, 1] % grid_size
    forces[:, 0] = Fx[gx, gy]
    forces[:, 1] = Fy[gx, gy]
    return forces

# -----------------------------
# Wrapper: compute forces (neighbor-list + gravity) with safety checks
# -----------------------------
def compute_forces(positions, velocities):
    Nloc = positions.shape[0]

    # Check for NaNs/Infs before building neighbor list
    if not cp.all(cp.isfinite(positions)):
        raise RuntimeError("Non-finite positions detected before neighbor build.")
    if not cp.all(cp.isfinite(velocities)):
        raise RuntimeError("Non-finite velocities detected before neighbor build.")

    # prepare cell arrays
    ncell = num_cells * num_cells
    cell_particles = cp.full((ncell * max_ppc,), -1, dtype=cp.int32)
    cell_counts = cp.zeros((ncell,), dtype=cp.int32)

    # launch neighbor build
    blocks = (Nloc + threads_per_block - 1) // threads_per_block
    try:
        build_neighbor_list_kernel((blocks,), (threads_per_block,),
                                  (positions.ravel(), cp.int32(Nloc),
                                   cp.float32(box_size), cp.float32(cell_size),
                                   cp.int32(num_cells),
                                   cell_particles, cell_counts, cp.int32(max_ppc)))
    except cp.cuda.runtime.CUDARuntimeError as e:
        raise RuntimeError(f"CUDA error building neighbor list: {e}")

    # prepare forces buffer (zeroed)
    forces = cp.zeros((Nloc, 2), dtype=cp.float32)

    # zero the flattened buffer used by kernel
    forces_flat = forces.ravel()

    cutoff2 = float(cell_size * cell_size)
    try:
        compute_local_forces_kernel((blocks,), (threads_per_block,),
                                    (positions.ravel(), velocities.ravel(), forces_flat,
                                     cp.int32(Nloc), cp.float32(box_size), cp.float32(cell_size),
                                     cp.int32(num_cells),
                                     cell_particles, cell_counts,
                                     cp.int32(max_ppc), cp.float32(epsilon), cp.float32(sigma),
                                     cp.float32(0.01), cp.float32(cutoff2), cp.float32(eps_safe)))
    except cp.cuda.runtime.CUDARuntimeError as e:
        raise RuntimeError(f"CUDA error during local force compute: {e}")

    # reshape forces_flat back to (N,2) view (forces already is that view)
    # add long-range gravity
    try:
        grav = compute_gravity_forces(positions)
    except Exception as e:
        raise RuntimeError(f"Error computing gravity forces: {e}")

    forces += grav
    return forces

# -----------------------------
# Integrator with safety clamps
# -----------------------------
def integrate(positions, velocities, forces):
    # standard velocity-Verlet step (GPU)
    positions = (positions + velocities * dt + 0.5 * forces * (dt * dt)) % box_size

    # ensure finite
    positions = cp.nan_to_num(positions, nan=0.0, posinf=0.0, neginf=0.0)

    # cap velocities to avoid runaways (after computing new forces)
    new_forces = compute_forces(positions, velocities)

    velocities = velocities + 0.5 * (forces + new_forces) * dt

    # clamp any NaN/Inf velocities to zero
    velocities = cp.nan_to_num(velocities, nan=0.0, posinf=v_max, neginf=-v_max)

    # velocity magnitude cap
    vnorm = cp.linalg.norm(velocities, axis=1, keepdims=True)
    too_fast = (vnorm > v_max).ravel()
    if cp.any(too_fast):
        # scale down only those exceeding v_max
        scale = (v_max / (vnorm + 1e-20))
        scale = cp.minimum(scale, 1.0)
        velocities = velocities * scale

    return positions, velocities, new_forces

# -----------------------------
# Save snapshot
# -----------------------------
def save_snapshot(positions, velocities, step, fname):
    try:
        cp.savez(fname, positions=positions, velocities=velocities, step=np.int64(step))
        print(f"Saved snapshot: {fname}")
    except Exception as e:
        print(f"Failed to save snapshot {fname}: {e}")

# -----------------------------
# Pygame visualization and main loop
# -----------------------------
def main_loop():
    pygame.init()
    screen = pygame.display.set_mode((600, 600))
    pygame.display.set_caption("FluidGravity (safe)")
    clock = pygame.time.Clock()

    # initial forces
    try:
        forces = compute_forces(positions, velocities)
    except Exception as e:
        print("Initial compute_forces failed:", e)
        save_snapshot(positions, velocities, -1, "crash_snapshot_initial.npz")
        raise

    for step in range(steps):
        try:
            for event in pygame.event.get():
                if event.type == QUIT:
                    save_snapshot(positions, velocities, step, "exit_snapshot.npz")
                    pygame.quit()
                    return

            # integrator (GPU)
            positions_local, velocities_local, forces_local = integrate(positions, velocities, forces)

            # update references (assign back to global arrays)
            # Note: we intentionally overwrite variables to keep memory contiguous
            positions[:] = positions_local
            velocities[:] = velocities_local
            forces[:] = forces_local

            # check for NaN/Infs and abort gracefully if detected
            if not cp.all(cp.isfinite(positions)) or not cp.all(cp.isfinite(velocities)):
                print("Non-finite values detected at step", step)
                save_snapshot(positions, velocities, step, "crash_snapshot_nonfinite.npz")
                raise RuntimeError("Non-finite values detected; aborting.")

            # render occasionally (reduce cost)
            if (step % 5) == 0:
                screen.fill((0, 0, 0))
                pos_np = cp.asnumpy(positions / box_size * 600.0).astype(int)
                for x, y in pos_np:
                    # clamp pixel coords inside screen
                    x = max(0, min(599, int(x)))
                    y = max(0, min(599, int(y)))
                    pygame.draw.circle(screen, (255, 255, 255), (x, y), 2)
                pygame.display.flip()
                clock.tick(60)

            # periodic saving
            if step % save_interval == 0 and step > 0:
                save_snapshot(positions, velocities, step, f"snapshot_step_{step}.npz")

        except Exception as e:
            # catch runtime/cuda errors, save state, and print traceback
            print("Exception at step", step, ":", e)
            traceback.print_exc()
            save_snapshot(positions, velocities, step, f"crash_snapshot_step_{step}.npz")
            pygame.quit()
            return

    # final save
    save_snapshot(positions, velocities, steps, f"final_snapshot_{steps}.npz")
    pygame.quit()

if __name__ == "__main__":
    try:
        main_loop()
    except Exception as e:
        print("Fatal error in main:", e)
        traceback.print_exc()
        # on fatal error try to save last-known state if available
        try:
            save_snapshot(positions, velocities, -999, "fatal_snapshot.npz")
        except Exception:
            pass
        sys.exit(1)
