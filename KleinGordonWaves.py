"""
kg2d_pygame.py

2D Klein-Gordon simulator with pygame visualization.

Equation: phi_tt = Laplacian(phi) - m^2 * phi  (units c = 1)
Numerical scheme: leap-frog (phi^{n+1} = 2 phi^n - phi^{n-1} + dt^2*(Lap phi^n - m^2 phi^n) )
Laplacian computed spectrally (FFT), so periodic BCs and easy mode-space editing.

Plug your mode-rounding code into apply_mode_rounding(phi_k).
"""

import numpy as np
from numpy.fft import fft2, ifft2, fftfreq, fftshift
import pygame
import time
import os

# ---------------------------
# Simulation parameters
# ---------------------------
N = 256                # grid points per side (power of two recommended)
L = 200               # physical box length
dx = L / N
m = 1                 # mass term (m^2 in KG)
dt = 0.1 * dx          # time step (CFL-like; reduce if unstable)
damping = 0          # optional small damping gamma: phi_tt + gamma phi_t = ...
t_max = 1e9            # dummy cap if wanted

# Visualization params
vis_scale = 1.0        # multiply phi to make contrast; tweak as needed
cmap = 'RdYlBu'        # not used directly; we build a simple colormap below
window_size = 768      # pygame window size in pixels (square)

# Mode rounding toggle (you will implement rounding in apply_mode_rounding)
mode_rounding_enabled = True

# Save snapshots directory
SNAPSHOT_DIR = "kg_snapshots"
os.makedirs(SNAPSHOT_DIR, exist_ok=True)

# ---------------------------
# Derived quantities & precompute
# ---------------------------
x = np.linspace(0, L, N, endpoint=False)
kx = 2 * np.pi * fftfreq(N, d=dx)   # frequency coordinates
ky = kx.copy()
KX, KY = np.meshgrid(kx, ky, indexing='ij')
k2 = KX**2 + KY**2

# Spectral operator for Laplacian: Lap phi -> ifft( -k^2 * phi_k )
L_spec = -k2

# Precompute mass term effect in update formula if desired
# (we compute in real space multiplied by dt^2)

# --------------------------------
# Helper: user-mode rounding hook
# --------------------------------
def apply_mode_rounding(phi_k):
    """
    Robust stochastic mode rounding.

    Each mode's |phi_k|^2 is treated as expected occupation number.
    For small amplitudes we Poisson-sample it.
    For large amplitudes (lam > 1e6), we approximate Poisson by a
    Gaussian with mean=lam and var=lam for numerical stability.
    """

    rho = np.abs(phi_k)**2
    mask = rho > 0.0

    rho0 = np.zeros_like(rho)

    # safe Poisson for small means
    small_mask = (mask) & (rho < 1e6)
    if np.any(small_mask):
        rho0[small_mask] = np.random.poisson(rho[small_mask])

    # Gaussian approximation for large lambda
    large_mask = (mask) & (~small_mask)
    if np.any(large_mask):
        lam = rho[large_mask]
        rho0[large_mask] = np.maximum(
            np.random.normal(lam, np.sqrt(lam)), 0.0
        )

    # reconstruct new complex amplitudes
    phi_k_new = np.zeros_like(phi_k, dtype=np.complex128)
    phi_k_phase = np.zeros_like(phi_k, dtype=np.complex128)
    phi_k_phase[mask] = phi_k[mask] / np.sqrt(rho[mask])
    phi_k_new[mask] = phi_k_phase[mask] * np.sqrt(rho0[mask])

    return phi_k_new


# --------------------------------
# Initialize fields
# --------------------------------
def init_fields():
    #rng = np.random.RandomState(seed)
    #phi = 1e-3 * rng.normal(size=(N, N))      # small random initial displacement
    #phi_prev = phi.copy()                     # start with zero velocity: phi_prev = phi - dt * phi_dot(0)
    # optionally give a localized bump:
    cx, cy = N//2, N//2
    r2 = (np.arange(N)[:,None]-cx)**2 + (np.arange(N)[None,:]-cy)**2
    phi = 0.5*np.exp(-r2/(2*(N*0.05)**2))
    phi_prev = np.zeros_like(phi)
    return phi, phi_prev

phi, phi_prev = init_fields()
phi_dot = np.zeros_like(phi)   # could be used if you want to include velocity explicitly

# ---------------------------
# Pygame color map helper
# ---------------------------
def build_colormap(n=256):
    # Build a simple diverging colormap similar to RdYlBu
    import colorsys
    cmap = np.zeros((n, 3), dtype=np.uint8)
    for i in range(n):
        t = i / (n-1)
        # map t in [0,1] to hue range 0.0 (blue) -> 0.0 (red) via a simple palette
        # We'll create blue -> white -> red style:
        if t < 0.5:
            # blue to white
            u = t / 0.5
            r = int(255 * (0.9 * u + 0.1 * (1-u)))
            g = int(255 * (0.9 * u + 0.1 * (1-u)))
            b = int(255 * (1.0 * (1-u) + 1.0 * u))
        else:
            # white to red
            u = (t-0.5) / 0.5
            r = int(255 * (1.0 * (1-u) + 1.0 * u))
            g = int(255 * (0.9 * (1-u) + 0.1 * u))
            b = int(255 * (0.9 * (1-u) + 0.1 * u))
        cmap[i] = np.clip([r,g,b], 0, 255)
    return cmap

CMAP = build_colormap(256)

def field_to_surface(phi_field, vmin=None, vmax=None):
    """
    Convert phi_field (2D array) to a pygame Surface (RGB).
    vmin/vmax control contrast. If None, autoscale using robust percentiles.
    """
    if vmin is None or vmax is None:
        p1, p99 = np.percentile(phi_field, [1, 99])
        vmin = p1 if vmin is None else vmin
        vmax = p99 if vmax is None else vmax
        if vmax - vmin < 1e-8: vmax = vmin + 1e-8
    norm = np.clip((phi_field - vmin) / (vmax - vmin), 0.0, 1.0)
    idx = (norm * 255).astype(np.int32)
    rgb = CMAP[idx]   # shape (N,N,3)
    # convert to pygame surface sized to window_size
    surf = pygame.surfarray.make_surface(np.flipud(rgb))  # flip so y orientation looks right
    surf = pygame.transform.smoothscale(surf, (window_size, window_size))
    return surf

# ---------------------------
# Simulation core: leap-frog step
# ---------------------------
def step_leapfrog(phi, phi_prev, dt, apply_rounding=True):
    """
    Perform a single leap-frog update returning phi_next, phi (new prev)
    phi_next = 2 phi - phi_prev + dt^2 * ( Laplacian(phi) - m^2 * phi ) - gamma * (phi - phi_prev)
    """
    # Spectral transform of phi
    phi_k = fft2(phi)

    # If user rounding is enabled, call hook (work in centered or raw ordering as you like)
    if mode_rounding_enabled and apply_rounding:
        phi_k = apply_mode_rounding(phi_k)

    # Compute Laplacian in spectral space: Lap phi -> ifft( -k^2 * phi_k )
    lap_phi = np.real(ifft2(L_spec * phi_k))

    # Leap-frog update
    phi_next = 2.0 * phi - phi_prev + dt*dt * (lap_phi - (m**2) * phi)

    return phi_next

# ---------------------------
# Init pygame
# ---------------------------
pygame.init()
screen = pygame.display.set_mode((window_size, window_size))
pygame.display.set_caption("2D Klein-Gordon (phi) — press SPACE to pause, r reset, c clear, s save snapshot")
clock = pygame.time.Clock()

paused = False
running = True
step_count = 0
last_snapshot = None

# Diagnostic: energy (approx)
def compute_energy(phi, phi_prev):
    # kinetic ~ ( (phi - phi_prev) / dt )^2 / 2
    phi_t = (phi - phi_prev) / dt
    kinetic = 0.5 * np.sum(phi_t**2) * (dx**2)
    # potential: 0.5*(|grad phi|^2 + m^2 phi^2)
    # compute grad in spectral space
    phi_k = fft2(phi)
    grad_x = np.real(ifft2(1j * KX * phi_k))
    grad_y = np.real(ifft2(1j * KY * phi_k))
    grad2 = grad_x**2 + grad_y**2
    potential = 0.5 * np.sum(grad2 + (m**2) * phi**2) * (dx**2)
    return kinetic + potential

# ---------------------------
# Main loop
# ---------------------------
try:
    t = 0.0
    t0 = time.time()
    while running:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                running = False
            elif ev.type == pygame.KEYDOWN:
                if ev.key == pygame.K_SPACE:
                    paused = not paused
                elif ev.key == pygame.K_r:
                    phi, phi_prev = init_fields()
                elif ev.key == pygame.K_c:
                    phi[:] = 0.0
                    phi_prev[:] = 0.0
                elif ev.key == pygame.K_s:
                    # save snapshot
                    surf = field_to_surface(phi * vis_scale)
                    fname = os.path.join(SNAPSHOT_DIR, f"phi_snapshot_{step_count:06d}.png")
                    pygame.image.save(surf, fname)
                    last_snapshot = fname
                    print("Saved snapshot:", fname)
                elif ev.key == pygame.K_q:
                    running = False
                elif ev.key == pygame.K_m:
                    # toggle mode rounding on/off
                    mode_rounding_enabled = not mode_rounding_enabled
                    print("Mode rounding:", mode_rounding_enabled)

        if not paused:
            # update
            phi_next = step_leapfrog(phi, phi_prev, dt, apply_rounding=True)
            phi_prev, phi = phi, phi_next
            t += dt
            step_count += 1

        # render
        surf = field_to_surface(phi * vis_scale)
        screen.blit(surf, (0,0))

        # overlay some info text
        font = pygame.font.SysFont("consolas", 16)
        txt = f"step {step_count} t={t:.3e} paused={paused} mode_round={mode_rounding_enabled}"
        energy = compute_energy(phi, phi_prev)
        txt2 = f"E≈{energy:.3e} m={m} dt={dt:.3e} dx={dx:.3e}"
        txt_surf = font.render(txt, True, (255,255,255))
        txt2_surf = font.render(txt2, True, (255,255,255))
        screen.blit(txt_surf, (6,6))
        screen.blit(txt2_surf, (6, 26))

        pygame.display.flip()
        clock.tick(60)   # limit to 60 FPS for UI; simulation may be faster internally

except KeyboardInterrupt:
    print("Interrupted by user")

finally:
    pygame.quit()
    if last_snapshot:
        print("Last snapshot:", last_snapshot)
    print("Exiting.")
