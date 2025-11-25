import numpy as np
import pygame
screen_size = 500
x0 = 5
y0 = -2
dr = 35
px0 = 0.1
py0 = 0.3
N = 500
M = 0.25
R = np.linspace(-N/2, N/2, N, endpoint=False)
P = np.linspace(-1/2, 1/2, N, endpoint=False)
X, Y= np.meshgrid(R, R, indexing='ij')
Px, Py= np.meshgrid(P, P, indexing='ij')
k_squared = Px**2 + Py**2
omega = np.sqrt(M**2 + k_squared)
K = np.sqrt(k_squared)
Psi = np.exp(-((X-x0)**2 + (Y-y0)**2)/(2*dr**2))*np.exp(-1j*(px0*X+py0*Y))/15
# Initialize Pygame
pygame.init()
screen = pygame.display.set_mode((screen_size, screen_size))
clock = pygame.time.Clock()
eps = 1/10**10
particles = []
Nparticles = np.zeros((N,N))
#def normalize_density(density):
#    """Normalize density for visualization"""
#    density = np.real(density)
#    return ((density+np.min(density)) / (np.max(density) + np.min(density))) * 255
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_SPACE:
                # Bosonic Field Mode Collapse Logic(Square Root; Poisson Distribution)
                PsiHat = np.fft.fft2(Psi)
                density1 = np.abs(PsiHat) ** 2
                NumPoints1 = np.random.poisson(density1)
                PsiHat = PsiHat * np.sqrt(NumPoints1/(density1 + eps))
                Psi = np.fft.ifft2(PsiHat)
                # Same Collapse Logic for Momentum/Position Domain
                density = np.abs(Psi) ** 2
                NumPoints = np.random.poisson(density)
                Psi = Psi * np.sqrt(NumPoints / (density + eps))
                # Particles
                Nparticles = NumPoints
                # Pilot Wave Trajectories Initialization
                particles = []
                for x in range(N):
                    for y in range(N):
                        count = int(Nparticles[x, y])
                        for _ in range(count):
                            particles.append([x, y])  # initial grid position
            else:
                Nparticles = np.zeros((N,N))
    # Klein Gordon matter Loop(Fourier Pseudo Spectral)
    PsiHat = np.fft.fft2(Psi)
    PsiHat = PsiHat * np.exp(1j * omega)
    Psi = np.fft.ifft2(PsiHat)
    # Pilot Wave Trajectories
    phase = np.angle(Psi)
    px = np.gradient(phase, axis=0)
    py = np.gradient(phase, axis=1)
    p_squared = px**2 + py**2
    E = np.sqrt(M ** 2 + p_squared)
    vx = px/E
    vy = py/E
    vt = 1
    for i, (x, y) in enumerate(particles):
        dx = vx[int(x) % N, int(y) % N]
        dy = vy[int(x) % N, int(y) % N]
        particles[i][0] = (x + dx) % N
        particles[i][1] = (y + dy) % N
    # Visual for the Spawned Particles(per every Collapse)
    density0 = Nparticles #np.sum(density, axis=2)
    img = (density0).astype(np.uint8)
    # Convert to RGB grayscale (each pixel: (val, val, val))
    gray_rgb = np.stack([img] * 3, axis=-1)
    # Now make surface properly
    surf = pygame.surfarray.make_surface(gray_rgb.swapaxes(0, 1))
    surf = pygame.transform.scale(surf, (screen_size, screen_size))
    screen.blit(surf, (0, 0))
    # Visualization
    for x, y in particles:
        pygame.draw.circle(screen, (255, 0, 0), (int(y * screen_size / N), int(x * screen_size / N)), 2)
    pygame.display.flip()
    clock.tick(60)
pygame.quit()

