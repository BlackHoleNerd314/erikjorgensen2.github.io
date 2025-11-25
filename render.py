# render.py

#class MyClass:
#    def __init__(self, value):
#        self.value = value  # public attribute
#
#    def display(self):
#        return self.value  # public method

import pygame
import sys
import numpy as np
from player import KeyInput,StepMove,Camera
from blackhole import OrbitPair

particles = 24

def CameraView(dim,xyz):
    r0 = xyz[0]
    X = xyz[1]
    Y = xyz[2]
    Z = xyz[3]
    cond0 = ((X/Y>6)or(X/Y<(-6)))or((Z/Y>4)or(Z/Y<(-4)))
    if (Y <= 0)or(cond0):
        Xpix = 1
        Ypix = 1
        r = 0
    else:
        Xpix = 60*(6+X/Y)
        Ypix = 60*(4-Z/Y)
        r = r0*60/Y
    if dim==3:
        return r
    elif dim==1:
        return Xpix
    elif dim==2:
        return Ypix


def Render():
    # Initialize pygame
    pygame.init()

    data0 = np.array((0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0))

    # Set up the display
    screen = pygame.display.set_mode((720, 480))
    pygame.display.set_caption("Cosmos Simulator")

    clock = pygame.time.Clock()

    # Main loop
    running = True
    while running:
        # Event handling
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        # Fill the screen with black
        screen.fill((0, 0, 0))

        data1 = data0
        key0 = KeyInput()
        data0 = StepMove(data1, key0)
        ###print(data0)
        # List of points (each is an (x, y) tuple)
        #points = [(100, 100), (200, 200), (300, 300), (400, 150)]

        #xyz = [1,0,10,0]#Camera(data0, [3,0,100,0,1])
        for N in range(0,1000):
            xyz = Camera(data0, [N*2, 0, 0.5, 0, 1])
            pygame.draw.circle(screen, (255, 255, 255), (CameraView(1,xyz),CameraView(2,xyz)), CameraView(3,xyz))

        # Update the display
        pygame.display.flip()
        clock.tick(60)

    # Quit pygame and exit the program
    pygame.quit()
    sys.exit()


