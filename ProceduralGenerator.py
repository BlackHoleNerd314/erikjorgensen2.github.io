import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi

# Recursive function to create fractal Voronoi structure using boundary edges for galaxy clusters
def fractal_voronoi(points, depth, target_depth):
    if depth == target_depth:
        return points
    new_points = []
    vor = Voronoi(points)
    for ridge in vor.ridge_vertices:
        if -1 not in ridge:
            pt1 = vor.vertices[ridge[0]]
            pt2 = vor.vertices[ridge[1]]
            # Sample points along the ridge to form a filament
            steps = 4
            for i in range(1, steps):
                t = i / steps
                interp = pt1 * (1 - t) + pt2 * t
                new_points.append(interp)
    return fractal_voronoi(np.array(new_points), depth + 1, target_depth)

# User input for recursion depth
iteration_count = int(input("Enter the number of fractal iterations (e.g. 3 to 5): "))

# Initial grid seed points with higher resolution
x = np.linspace(0, 300, 10)
y = np.linspace(0, 300, 10)
xv, yv = np.meshgrid(x, y)
initial_points = np.vstack([xv.ravel(), yv.ravel()]).T

# Generate only the final iteration's points
final_points = fractal_voronoi(initial_points, depth=0, target_depth=iteration_count)

# Plot the most detailed cosmic web structure
plt.figure(figsize=(10, 10), facecolor='black')
plt.scatter(final_points[:, 0], final_points[:, 1], s=0.5, color='white')
plt.title("2D Cosmic Web Fractal (Most Detailed Iteration)", color='white')
plt.gca().set_aspect('equal', adjustable='box')
plt.axis('off')
plt.gca().set_facecolor('black')
plt.xlim(0, 300)
plt.ylim(0, 300)
plt.show()