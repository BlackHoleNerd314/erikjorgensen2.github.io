import math
import numpy as np
import pygame
import sys

xyz20 = [
390.0,	3.769647E-03,	4.146161E-04,	1.847260E-02,
395.0,	9.382967E-03,	1.059646E-03,	4.609784E-02,
400.0,	2.214302E-02,	2.452194E-03,	1.096090E-01,
405.0,	4.742986E-02,	4.971717E-03,	2.369246E-01,
410.0,	8.953803E-02,	9.079860E-03,	4.508369E-01,
415.0,	1.446214E-01,	1.429377E-02,	7.378822E-01,
420.0,	2.035729E-01,	2.027369E-02,	1.051821E+00,
425.0,	2.488523E-01,	2.612106E-02,	1.305008E+00,
430.0,	2.918246E-01,	3.319038E-02,	1.552826E+00,
435.0,	3.227087E-01,	4.157940E-02,	1.748280E+00,
440.0,	3.482554E-01,	5.033657E-02,	1.917479E+00,
445.0,	3.418483E-01,	5.743393E-02,	1.918437E+00,
450.0,	3.224637E-01,	6.472352E-02,	1.848545E+00,
455.0,	2.826646E-01,	7.238339E-02,	1.664439E+00,
460.0,	2.485254E-01,	8.514816E-02,	1.522157E+00,
465.0,	2.219781E-01,	1.060145E-01,	1.428440E+00,
470.0,	1.806905E-01,	1.298957E-01,	1.250610E+00,
475.0,	1.291920E-01,	1.535066E-01,	9.991789E-01,
480.0,	8.182895E-02,	1.788048E-01,	7.552379E-01,
485.0,	4.600865E-02,	2.064828E-01,	5.617313E-01,
490.0,	2.083981E-02,	2.379160E-01,	4.099313E-01,
495.0,	7.097731E-03,	2.850680E-01,	3.105939E-01,
500.0,	2.461588E-03,	3.483536E-01,	2.376753E-01,
505.0,	3.649178E-03,	4.277595E-01,	1.720018E-01,
510.0,	1.556989E-02,	5.204972E-01,	1.176796E-01,
515.0,	4.315171E-02,	6.206256E-01,	8.283548E-02,
520.0,	7.962917E-02,	7.180890E-01,	5.650407E-02,
525.0,	1.268468E-01,	7.946448E-01,	3.751912E-02,
530.0,	1.818026E-01,	8.575799E-01,	2.438164E-02,
535.0,	2.405015E-01,	9.071347E-01,	1.566174E-02,
540.0,	3.098117E-01,	9.544675E-01,	9.846470E-03,
545.0,	3.804244E-01,	9.814106E-01,	6.131421E-03,
550.0,	4.494206E-01,	9.890228E-01,	3.790291E-03,
555.0,	5.280233E-01,	9.994608E-01,	2.327186E-03,
560.0,	6.133784E-01,	9.967737E-01,	1.432128E-03,
565.0,	7.016774E-01,	9.902549E-01,	8.822531E-04,
570.0,	7.967750E-01,	9.732611E-01,	5.452416E-04,
575.0,	8.853376E-01,	9.424569E-01,	3.386739E-04,
580.0,	9.638388E-01,	8.963613E-01,	2.117772E-04,
585.0,	1.051011E+00,	8.587203E-01,	1.335031E-04,
590.0,	1.109767E+00,	8.115868E-01,	8.494468E-05,
595.0,	1.143620E+00,	7.544785E-01,	5.460706E-05,
600.0,	1.151033E+00,	6.918553E-01,	3.549661E-05,
605.0,	1.134757E+00,	6.270066E-01,	2.334738E-05,
610.0,	1.083928E+00,	5.583746E-01,	1.554631E-05,
615.0,	1.007344E+00,	4.895950E-01,	1.048387E-05,
620.0,	9.142877E-01,	4.229897E-01,	0.000000E+00,
625.0,	8.135565E-01,	3.609245E-01,	0.000000E+00,
630.0,	6.924717E-01,	2.980865E-01,	0.000000E+00,
635.0,	5.755410E-01,	2.416902E-01,	0.000000E+00,
640.0,	4.731224E-01,	1.943124E-01,	0.000000E+00,
645.0,	3.844986E-01,	1.547397E-01,	0.000000E+00,
650.0,	2.997374E-01,	1.193120E-01,	0.000000E+00,
655.0,	2.277792E-01,	8.979594E-02,	0.000000E+00,
660.0,	1.707914E-01,	6.671045E-02,	0.000000E+00,
665.0,	1.263808E-01,	4.899699E-02,	0.000000E+00,
670.0,	9.224597E-02,	3.559982E-02,	0.000000E+00,
675.0,	6.639960E-02,	2.554223E-02,	0.000000E+00,
680.0,	4.710606E-02,	1.807939E-02,	0.000000E+00,
685.0,	3.292138E-02,	1.261573E-02,	0.000000E+00,
690.0,	2.262306E-02,	8.661284E-03,	0.000000E+00,
695.0,	1.575417E-02,	6.027677E-03,	0.000000E+00,
700.0,	1.096778E-02,	4.195941E-03,	0.000000E+00,
705.0,	7.608750E-03,	2.910864E-03,	0.000000E+00,
710.0,	5.214608E-03,	1.995557E-03,	0.000000E+00,
715.0,	3.569452E-03,	1.367022E-03,	0.000000E+00,
720.0,	2.464821E-03,	9.447269E-04,	0.000000E+00,
725.0,	1.703876E-03,	6.537050E-04,	0.000000E+00,
730.0,	1.186238E-03,	4.555970E-04,	0.000000E+00,
735.0,	8.269535E-04,	3.179738E-04,	0.000000E+00,
740.0,	5.758303E-04,	2.217445E-04,	0.000000E+00,
745.0,	4.058303E-04,	1.565566E-04,	0.000000E+00,
750.0,	2.856577E-04,	1.103928E-04,	0.000000E+00,
755.0,	2.021853E-04,	7.827442E-05,	0.000000E+00,
760.0,	1.438270E-04,	5.578862E-05,	0.000000E+00,
765.0,	1.024685E-04,	3.981884E-05,	0.000000E+00,
770.0,	7.347551E-05,	2.860175E-05,	0.000000E+00,
775.0,	5.259870E-05,	2.051259E-05,	0.000000E+00,
780.0,	3.806114E-05,	1.487243E-05,	0.000000E+00,
785.0,	2.758222E-05,	1.080001E-05,	0.000000E+00,
790.0,	2.004122E-05,	7.863920E-06,	0.000000E+00,
795.0,	1.458792E-05,	5.736935E-06,	0.000000E+00,
800.0,	1.068141E-05,	4.211597E-06,	0.000000E+00,
805.0,	7.857521E-06,	3.106561E-06,	0.000000E+00,
810.0,	5.768284E-06,	2.286786E-06,	0.000000E+00,
815.0,	4.259166E-06,	1.693147E-06,	0.000000E+00,
820.0,	3.167765E-06,	1.262556E-06,	0.000000E+00,
825.0,	2.358723E-06,	9.422514E-07,	0.000000E+00,
830.0,	1.762465E-06,	7.053860E-07,	0.000000E+00]
x = np.zeros((1000,1))
y = np.zeros((1000,1))
z = np.zeros((1000,1))
xyz1 = np.zeros((88,4))
for ind0 in range(0,88):
    xyz1[ind0,0] = xyz20[ind0*4+0]
    xyz1[ind0, 1] = xyz20[ind0*4 + 1]
    xyz1[ind0, 2] = xyz20[ind0 * 4 + 2]
    xyz1[ind0, 3] = xyz20[ind0 * 4 + 3]
for Vind in range(78,166):
    Vind0 = Vind - 78
    x[Vind] = xyz1[Vind0,1]
    y[Vind] = xyz1[Vind0,2]
    z[Vind] = xyz1[Vind0,3]
for UVind in range(0,78):
    UVind0 = 78-UVind
    Yfact = y[78]/y[78+1]
    x[UVind0] = x[UVind0+1]*Yfact
    y[UVind0] = y[UVind0+1]*Yfact
    z[UVind0] = z[UVind0+1]*Yfact
for Rind in range(123,166):
    Zfact = z[123]/z[123-1]
    z[Rind] = z[Rind-1]*Zfact
for IRind in range(166,1000-1):
    Yfact0 = y[165]/y[164]
    x[IRind] = x[IRind-1]*Yfact0
    y[IRind] = y[IRind-1]*Yfact0
    z[IRind] = z[IRind-1]*Yfact0
fix = 1/10**24
xy0 = np.zeros((1000,1))
for ind0 in range(0,1000):
    xy0[ind0] = x[ind0]/(y[ind0]+fix)
xy = min(xy0)
for ind0 in range(0,1000):
    x[ind0] = ((xy0[ind0])-xy)*y[ind0]
x = x/sum(x)
y = y/sum(y)
z = z/sum(z)
r0 = x
b0 = z
g0 = (y-0.3826*x+0.0559*z)/(1-0.3826+0.0559)
g = g0/sum(g0)
r = r0/sum(r0)
b = b0/sum(b0)
nm0 = np.linspace(0,1000,1000)
nm1 = 0
y0 = np.zeros((1000,1))
for iter0 in nm0:
    r0[nm1] = r[nm1]*(nm1**4)
    g0[nm1] = g[nm1]*(nm1**4)
    b0[nm1] = b[nm1]*(nm1**4)
    y0[nm1] = y[nm1]*(nm1**4)
    nm1 += 1
r = np.sqrt(r0)
g = np.sqrt(g0)
b = np.sqrt(b0)
y = np.sqrt(y0)
I = np.zeros((5,5,3))
RGBmax0 = max([max(r),max(g),max(b),max(y)])
R = np.round((r/RGBmax0)*255)
G = np.round((g/RGBmax0)*255)
B = np.round((b/RGBmax0)*255)
L = np.round((y/RGBmax0)*255)
#plt.plot(nm0)
#plt.show()
# Initialize Pygame
pygame.init()
import pygame
import numpy as np

# Initialize Pygame mixer
pygame.mixer.pre_init(frequency=44100, size=-16, channels=1)




# Set up display
width, height = 720, 480
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("The Visible Light Spectrum")
clock = pygame.time.Clock()
t0 = 1
fps = 25
N = 0
N0 = 1
frame = 1

# New state variables for audio control
current_N = None # Tracks the N value of the currently playing sound
current_sound = None # Holds the currently playing pygame.mixer.Sound object
buffer_duration = 1.0 # Generate a 1-second buffer for seamless looping


# Main loop
running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    frame = frame + 1
    keys = pygame.key.get_pressed()
    if keys[pygame.K_UP]:
        N0 = N0 * 1.1
    elif keys[pygame.K_DOWN]:
        N0 = N0 / 1.1
    pygame.display.flip()
    if keys[pygame.K_RIGHT]:
        N = N + N0
    elif keys[pygame.K_LEFT]:
        N = N - N0



    # Parameters
    frequency = 440 * 2 ** (N / 12)  # *(1/(4.35456*10**17))
    duration = 0.04  # seconds
    sample_rate = 44100  # Hz
    amplitude = 32767  # max for 16-bit audio
    # Generate time array
    t = frame * duration + np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    # Generate sine wave
    waveform = (amplitude * np.sin(2 * np.pi * frequency * t)).astype(np.int16)
    # Convert NumPy array to Pygame Sound
    waveform_stereo = np.column_stack((waveform, waveform))
    sound = pygame.sndarray.make_sound(waveform_stereo)
    # Play sound
    sound.play()




    wavelength0 = (299792458/frequency)/(5*10**-9)

    ###wavelength0 = np.exp(N/N0)
    #wavelength0 = wavelength0 + 1
    # Define colors
    wavelength = round(wavelength0)
    wavelength1 = wavelength0*5/10**9

    BLACK = (0, 0, 0)
    screen.fill(BLACK)
    if wavelength>999:
        COLOR = BLACK
        WHITE = BLACK
    else:
        COLOR = (R[wavelength], G[wavelength], B[wavelength])
        WHITE = (L[wavelength],L[wavelength],L[wavelength])
    clock.tick(fps)

    screen.fill(COLOR)
    #sound.play()
    #pygame.time.wait(int(1 + duration * 1000))
    font = pygame.font.SysFont('Comic Sans', 24)
    # Or use a specific font file:
    # font = pygame.font.Font('path_to_font.ttf', 48)
    #text_surface = font.render(f"Wavelength in nanometers: {round(wavelength0*5)}", True, COLOR)  # White text
    #screen.blit(text_surface, (0, 0))  # Position (x=100, y=100)
    pygame.display.flip()  # Or pygame.display.update()
# Quit Pygame
pygame.mixer.quit()
pygame.quit()
sys.exit()