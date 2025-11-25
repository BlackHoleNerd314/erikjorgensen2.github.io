# main.py
import sys
import pygame

WIDTH, HEIGHT = 720, 480
FPS = 12

def main():
    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT))
    pygame.display.set_caption("Science Engine")
    clock = pygame.time.Clock()

    running = True
    while running:
        # Event handling
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            # Allow quitting with Esc
            elif event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE:
                running = False

        # Clear screen (black)
        screen.fill((0, 0, 0))

        # (Nothing drawn — empty window)

        pygame.display.flip()
        clock.tick(FPS)

    pygame.quit()
    sys.exit()

if __name__ == "__main__":
    main()
