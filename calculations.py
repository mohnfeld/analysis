# calculations.py

import numpy as np
from shapely.geometry import Point, LineString


def calculate_speed(x1, z1, t1, x2, z2, t2):
    distance = np.sqrt((x2 - x1) ** 2 + (z2 - z1) ** 2)
    time_diff = t2 - t1
    if time_diff == 0:
        return 0
    return distance / time_diff


def calculate_area_covered_and_speeds(pos_x, pos_z, speeds, time, outer_polygon, inner_polygon):
    trajectory_points = [Point(x, z) for x, z in zip(pos_x, pos_z)]

    # Calculate area coverage
    trajectory_line = LineString(zip(pos_x, pos_z))
    trajectory_area = trajectory_line.buffer(0.1)

    total_intersection = trajectory_area.intersection(outer_polygon)
    green_intersection = trajectory_area.intersection(inner_polygon)
    red_intersection = total_intersection.difference(green_intersection)

    total_area = outer_polygon.area
    green_area = inner_polygon.area
    red_area = total_area - green_area

    total_percentage_covered = (total_intersection.area / total_area) * 100
    green_percentage_covered = (green_intersection.area / green_area) * 100
    red_percentage_covered = (red_intersection.area / red_area) * 100

    # Calculate speeds, accelerations, stops, and time spent for safe and unsafe areas
    safe_speeds = []
    unsafe_speeds = []
    safe_accelerations = []
    unsafe_accelerations = []
    time_in_safe_area = 0
    time_in_unsafe_area = 0

    stop_threshold = 0.01  # m/s
    stop_duration_threshold = 0.5  # seconds

    safe_stops = []
    unsafe_stops = []
    current_stop_start = None
    current_stop_area = None

    # Calculate shortest distances to walls
    outer_coords = list(outer_polygon.exterior.coords)
    distances_to_walls = []  # Placeholder if needed
    average_distance_to_walls = None  # Placeholder if needed

    for i, (point, speed) in enumerate(zip(trajectory_points, speeds)):
        if inner_polygon.contains(point):
            safe_speeds.append(speed)
            if i > 0:
                time_diff = time[i] - time[i - 1]
                time_in_safe_area += time_diff
                if i > 1:
                    acceleration = (speed - speeds[i - 1]) / time_diff
                    safe_accelerations.append(acceleration)

            if speed < stop_threshold:
                if current_stop_start is None:
                    current_stop_start = time[i]
                    current_stop_area = 'safe'
            else:
                if current_stop_start is not None and current_stop_area == 'safe':
                    stop_duration = time[i] - current_stop_start
                    if stop_duration >= stop_duration_threshold:
                        safe_stops.append(stop_duration)
                    current_stop_start = None
                    current_stop_area = None

        elif outer_polygon.contains(point):
            unsafe_speeds.append(speed)
            if i > 0:
                time_diff = time[i] - time[i - 1]
                time_in_unsafe_area += time_diff
                if i > 1:
                    acceleration = (speed - speeds[i - 1]) / time_diff
                    unsafe_accelerations.append(acceleration)

            if speed < stop_threshold:
                if current_stop_start is None:
                    current_stop_start = time[i]
                    current_stop_area = 'unsafe'
            else:
                if current_stop_start is not None and current_stop_area == 'unsafe':
                    stop_duration = time[i] - current_stop_start
                    if stop_duration >= stop_duration_threshold:
                        unsafe_stops.append(stop_duration)
                    current_stop_start = None
                    current_stop_area = None

    # Check for any ongoing stop at the end of the trajectory
    if current_stop_start is not None:
        stop_duration = time[-1] - current_stop_start
        if stop_duration >= stop_duration_threshold:
            if current_stop_area == 'safe':
                safe_stops.append(stop_duration)
            else:
                unsafe_stops.append(stop_duration)

    # Calculate average distance to walls if needed
    # distances_to_walls = [shortest_distance_to_walls((x, z), outer_coords) for x, z in zip(pos_x, pos_z)]
    # average_distance_to_walls = np.mean(distances_to_walls)

    return (total_percentage_covered, red_percentage_covered, green_percentage_covered,
            safe_speeds, unsafe_speeds, safe_accelerations, unsafe_accelerations,
            time_in_safe_area, time_in_unsafe_area, safe_stops, unsafe_stops)
