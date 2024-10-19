# plotting_utils.py

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import pandas as pd
import os
from geometry_utils import rotate_points, calculate_tilt_angle
from calculations import calculate_speed, calculate_area_covered_and_speeds
from geometry_utils import shortest_distance_to_walls

def create_path_from_polygon(polygon):
    if polygon.is_empty:
        return None

    vertices = []
    codes = []

    def extract_coords(poly):
        # Exterior ring
        x, y = poly.exterior.coords.xy
        coords = list(zip(x, y))
        vertices.extend(coords)
        codes.extend([Path.MOVETO] + [Path.LINETO]*(len(coords)-2) + [Path.CLOSEPOLY])

        # Interior rings (holes)
        for interior in poly.interiors:
            x, y = interior.coords.xy
            coords = list(zip(x, y))
            vertices.extend(coords)
            codes.extend([Path.MOVETO] + [Path.LINETO]*(len(coords)-2) + [Path.CLOSEPOLY])

    if isinstance(polygon, Polygon):
        extract_coords(polygon)
    else:
        for poly in polygon.geoms:
            extract_coords(poly)

    return Path(vertices, codes)

def draw_elevated_platform(ax, angle, center_x, center_z):
    # Platform coordinates
    rotated_outer_coords = [(-6.00, 1.73), (-6.60, 1.73), (-6.60, 2.89), (-7.20, 2.89),
                            (-7.20, 4.66), (-7.80, 4.66), (-7.80, 7.00), (-4.80, 7.00),
                            (-4.80, 4.66), (-5.40, 4.66), (-5.40, 2.89), (-6.00, 2.89)]
    rotated_outer_coords.append(rotated_outer_coords[0])  # Close the polygon
    xs, zs = zip(*rotated_outer_coords)
    ax.plot(xs, zs, 'r-', linewidth=2)

    rotated_inner_coords = [(-6.78, 3.31), (-6.78, 5.08), (-7.38, 5.08), (-7.38, 6.58),
                            (-5.22, 6.58), (-5.22, 5.08), (-5.82, 5.08), (-5.82, 3.31)]
    rotated_inner_coords.append(rotated_inner_coords[0])  # Close the polygon
    inner_xs, inner_zs = zip(*rotated_inner_coords)
    ax.plot(inner_xs, inner_zs, 'g-', linewidth=2)  # Use green color for inner platform

    # Create shapely Polygons
    outer_polygon = Polygon(rotated_outer_coords)
    inner_polygon = Polygon(rotated_inner_coords)

    # Compute the area between the outer and inner polygons
    red_area_polygon = outer_polygon.difference(inner_polygon)
    green_area_polygon = inner_polygon

    # Create and plot patches
    red_area_path = create_path_from_polygon(red_area_polygon)
    if red_area_path:
        red_patch = PathPatch(red_area_path, facecolor='red', edgecolor='none', alpha=0.3)
        ax.add_patch(red_patch)

    inner_path = create_path_from_polygon(inner_polygon)
    green_patch = PathPatch(inner_path, facecolor='green', edgecolor='none', alpha=0.3)
    ax.add_patch(green_patch)

    # Plot center and pillar
    center_x_rotated, center_z_rotated = rotate_points(-5.75, 4.06, -angle, center_x, center_z)
    ax.plot(center_x_rotated, center_z_rotated, 'x')

    # Pillar coordinates
    pillar_coords = [(-6.26, 5.14), (-6.075, 5.19), (-6.03, 5.02), (-6.215, 4.97)]
    rotated_pillar = [rotate_points(x, z, -angle, center_x, center_z) for x, z in pillar_coords]
    rotated_pillar.append(rotated_pillar[0])  # Close the polygon
    pxs, pzs = zip(*rotated_pillar)
    ax.plot(pxs, pzs, 'r-')

    return outer_polygon, inner_polygon

def draw_training_room(ax, angle, center_x, center_z):
    platform_coords = [(-4.79, 5.74), (-5.67, 9.17), (-8.89, 8.32), (-8.00, 4.91)]
    rotated_coords = [rotate_points(x, z, angle, center_x, center_z) for x, z in platform_coords]
    rotated_coords.append(rotated_coords[0])
    xs, zs = zip(*rotated_coords)
    ax.plot(xs, zs, 'r-', linewidth=2)

def draw_empty_room(ax):
    pass

def draw_dark_maze_room(ax):
    coordinates = {
        "pillar": [(-6.1, 5.23), (-6.04, 5.05), (-6.23, 5.02), (-6.28, 5.18)],
        "person1": [(-6.76, 2.4), (-6.67, 2.16), (-7.21, 2.04), (-7.27, 2.3)],
        "person2": [(-6.68, 4.76), (-6.55, 4.26), (-6.84, 4.2), (-6.94, 4.72)],
        "person3": [(-7.51, 5.48), (-7.45, 5.15), (-8.01, 5.02), (-8.07, 5.34)],
        "wall1": [(-4.99, 2.49), (-4.94, 2.32), (-7.28, 1.69), (-7.33, 1.86)],
        "wall2": [(-5.56, 3.22), (-5.74, 3.17), (-5.96, 3.98), (-5.04, 4.21),
                  (-5.64, 6.45), (-6.29, 6.29), (-6.35, 6.49), (-5.61, 6.69),
                  (-4.93, 4.05), (-5.7, 3.84)],
        "wall3": [(-5.08, 8.22), (-5.02, 8.03), (-7.56, 7.37), (-7.62, 7.56)],
        "wall4": [(-7.26, 5.73), (-7.45, 5.67), (-7.27, 4.9), (-8.02, 4.72),
                  (-7.93, 4.51), (-7.22, 4.72), (-7.03, 4.08), (-6.84, 4.15),
                  (-7.0, 4.77), (-6.02, 5.01), (-6.08, 5.22), (-7.06, 4.98)]
    }
    angle = np.radians(-14.9)
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle),  np.cos(angle)]])
    center_point = np.array([-6.22, 5.005])

    rotated_coords = {}
    for k, v in coordinates.items():
        new_points = []
        for p in v:
            coord = np.array(p)
            new_point = np.dot(rotation_matrix, (coord - center_point)) + center_point
            new_points.append(new_point)
        new_points.append(new_points[0])  # Close the polygon
        rotated_coords[k] = new_points
        xs, zs = zip(*new_points)
        ax.plot(xs, zs, 'r-', linewidth=2)
    return rotated_coords

def plot_movement_trajectory_with_layout(file_path):
    # Updated to handle relative paths and cross-platform compatibility

    # Ensure the file_path is relative to the current directory
    file_path = os.path.join(os.getcwd(), file_path)

    # Extract scene ID from file name
    file_name = os.path.basename(file_path)
    scene_id = file_name.split('_')[2].split('.')[0]

    # Read movement data
    movement_data = pd.read_csv(file_path)
    pos_x = movement_data['pos_x']
    pos_z = movement_data['pos_z']
    time = movement_data['time']

    # Calculate tilt angle
    platform_coords = [(-4.83, 1.89), (-5.46, 1.71)]
    tilt_angle = calculate_tilt_angle(platform_coords)
    print(f"Calculated tilt angle: {np.degrees(tilt_angle):.2f} degrees")

    # Calculate center for rotation
    center_x = np.mean(pos_x)
    center_z = np.mean(pos_z)

    # Rotate positions
    pos_x_rotated, pos_z_rotated = rotate_points(pos_x, pos_z, -tilt_angle, center_x, center_z)

    # Calculate speeds
    speeds = [0]  # First point has no speed
    for i in range(1, len(pos_x_rotated)):
        speed = calculate_speed(pos_x_rotated[i-1], pos_z_rotated[i-1], time[i-1],
                                pos_x_rotated[i], pos_z_rotated[i], time[i])
        speeds.append(speed)

    # Plotting
    fig, ax1 = plt.subplots(figsize=(10, 10))
    sc = ax1.scatter(pos_x_rotated, pos_z_rotated, c=time, cmap='viridis', s=5)
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Z Position')
    ax1.set_title(f'Top-Down View of Straightened Movement Trajectory for {scene_id}')
    plt.colorbar(sc, ax=ax1, label='Time')
    ax1.grid(False)
    ax1.axis('equal')

    # Draw the appropriate layout based on the scene ID
    if scene_id == 'T011':
        outer_polygon, inner_polygon = draw_elevated_platform(ax1, tilt_angle, center_x, center_z)

        # Calculate area coverage, speeds, accelerations, stops, and time spent for safe/unsafe areas
        (total_area_covered, red_area_covered, green_area_covered,
         safe_speeds, unsafe_speeds, safe_accelerations, unsafe_accelerations,
         time_in_safe_area, time_in_unsafe_area, safe_stops, unsafe_stops,
         average_distance_to_walls) = calculate_area_covered_and_speeds(
            pos_x_rotated, pos_z_rotated, speeds, time, outer_polygon, inner_polygon)

        print(f"Total platform area covered: {total_area_covered:.2f}%")
        print(f"Unsafe (red) area covered: {red_area_covered:.2f}%")
        print(f"Safe (green) area covered: {green_area_covered:.2f}%")

        print(f"\nTime spent in safe area: {time_in_safe_area:.2f} seconds")
        print(f"Time spent in unsafe area: {time_in_unsafe_area:.2f} seconds")

        print(f"\nOverall speed statistics:")
        print(f"Min speed: {min(speeds):.2f}")
        print(f"Max speed: {max(speeds):.2f}")
        print(f"Avg speed: {np.mean(speeds):.2f}")

        # Additional statistics can be printed as needed

    elif scene_id == 'T003':
        draw_training_room(ax1, tilt_angle, center_x, center_z)
    elif scene_id == 'T007':
        draw_empty_room(ax1)
    elif scene_id == 'T015':
        coords = draw_dark_maze_room(ax1)
        # Optionally calculate mean minimum corner distance

    plt.tight_layout()
    plt.show()
