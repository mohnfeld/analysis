import os
import csv
import numpy as np
import pandas as pd

def calculate_tilt_angle(platform_coords):
    x1, z1 = platform_coords[0]
    x2, z2 = platform_coords[1]
    dx = x2 - x1
    dz = z2 - z1
    slope = dz / dx
    angle = np.arctan(slope)
    angle_deg = np.degrees(angle)
    return np.radians(-angle_deg)

def rotate_points(x, z, angle, center_x, center_z):
    x_centered = x - center_x
    z_centered = z - center_z
    x_rotated = x_centered * np.cos(angle) + z_centered * np.sin(angle)
    z_rotated = -x_centered * np.sin(angle) + z_centered * np.cos(angle)
    x_final = x_rotated + center_x
    z_final = z_rotated + center_z
    return x_final, z_final

if __name__ == '__main__':

    platform_coords =  platform_coords = [
        (-4.83, 1.89), (-5.46, 1.71), (-5.74, 2.84), (-6.32, 2.68),
        (-6.81, 4.40), (-7.39, 4.24), (-8.05, 6.47), (-5.11, 7.31),
        (-4.45, 5.06), (-5.03, 4.90), (-4.56, 3.17), (-5.14, 3.01)
    ]
    maze_coords = [(-7.28, 1.69), (-4.94, 2.32)]
    angles = {}
    angles["platform"] = calculate_tilt_angle(platform_coords)
    angles["maze"] = calculate_tilt_angle(maze_coords)
    print([np.degrees(a) for a in angles.values()])

    srcdir = r"C:\Users\lal\Documents\tez\analysis\BEA_edges_11.09.24\S001\trackers"
    targetdir = "C:\\Users\\lal\\Documents\\tez\\analysis\\BEA_edges_11.09.24\\S001\\trackers_rotated\\"

    for scene in ["platform", "maze"]:
        for f in os.listdir(os.path.join(srcdir,scene)):
            srcfile = os.path.abspath(os.path.join(srcdir, scene, f))
            targetfile = os.path.abspath(os.path.join(targetdir, scene, f))
            movement_data = pd.read_csv(srcfile)
            pos_x = movement_data['pos_x']
            pos_z = movement_data['pos_z']
            time = movement_data['time']
            center_x = np.mean(pos_x)
            center_z = np.mean(pos_z)
            movement_data["pos_x"], movement_data["pos_z"] = rotate_points(movement_data["pos_x"], movement_data["pos_z"], angles[scene], center_x, center_z)
            movement_data.to_csv(targetfile)
                