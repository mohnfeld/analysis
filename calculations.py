# calculations.py

import numpy as np
from shapely.geometry import Point, LineString


def df_diff(df):
    diff = df.diff(axis=0).abs()
    diff[diff > 180] = 360-diff[diff > 180]
    return diff
def calculate_total_rotation(rot_x, rot_y, rot_z):
    #
    angles = list(zip(*(rot_x, rot_y, rot_z)))
    angle_diff = np.abs(np.diff(angles,axis=0))
    angle_diff[angle_diff > 180] = 360-angle_diff[angle_diff > 180]
    total_rotation = np.mean(angle_diff,axis=0)
    return total_rotation

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
    total_distance= trajectory_line.length
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

#elevated platform functions


def calculate_trajectory_coverage(pos_x, pos_z, outer_polygon, inner_polygon):
    """Calculate coverage percentages for different areas for elevated platform."""
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

    return total_percentage_covered, red_percentage_covered, green_percentage_covered

def calculate_speeds_by_area(trajectory_points, speeds, inner_polygon, outer_polygon):
    """Calculate speeds in safe and unsafe areas."""
    safe_speeds = []
    unsafe_speeds = []
    
    for point, speed in zip(trajectory_points, speeds):
        if inner_polygon.contains(point):
            safe_speeds.append(speed)
        elif outer_polygon.contains(point):
            unsafe_speeds.append(speed)
            
    return safe_speeds, unsafe_speeds

def calculate_accelerations(speeds, time, trajectory_points, inner_polygon, outer_polygon):
    """Calculate accelerations in safe and unsafe areas."""
    safe_accelerations = []
    unsafe_accelerations = []
    
    for i in range(1, len(speeds)):
        time_diff = time[i] - time[i-1]
        acceleration = (speeds[i] - speeds[i-1]) / time_diff
        point = trajectory_points[i]
        
        if inner_polygon.contains(point):
            safe_accelerations.append(acceleration)
        elif outer_polygon.contains(point):
            unsafe_accelerations.append(acceleration)
            
    return safe_accelerations, unsafe_accelerations

def calculate_time_in_areas(trajectory_points, time, inner_polygon, outer_polygon):
    """Calculate time spent in safe and unsafe areas."""
    time_in_safe_area = 0
    time_in_unsafe_area = 0
    
    for i in range(1, len(trajectory_points)):
        time_diff = time[i] - time[i-1]
        point = trajectory_points[i]
        
        if inner_polygon.contains(point):
            time_in_safe_area += time_diff
        elif outer_polygon.contains(point):
            time_in_unsafe_area += time_diff
            
    return time_in_safe_area, time_in_unsafe_area

def calculate_stops(trajectory_points, speeds, time, inner_polygon, outer_polygon, 
                   stop_threshold=0.01, stop_duration_threshold=0.5):
    """Calculate stops in safe and unsafe areas."""
    safe_stops = []
    unsafe_stops = []
    current_stop_start = None
    current_stop_area = None

    for i, (point, speed) in enumerate(zip(trajectory_points, speeds)):
        is_safe = inner_polygon.contains(point)
        is_unsafe = outer_polygon.contains(point) if not is_safe else False
        
        if speed < stop_threshold:
            if current_stop_start is None:
                current_stop_start = time[i]
                current_stop_area = 'safe' if is_safe else 'unsafe'
        else:
            if current_stop_start is not None:
                stop_duration = time[i] - current_stop_start
                if stop_duration >= stop_duration_threshold:
                    if current_stop_area == 'safe':
                        safe_stops.append(stop_duration)
                    else:
                        unsafe_stops.append(stop_duration)
                current_stop_start = None
                current_stop_area = None

    # Check for ongoing stop at the end
    if current_stop_start is not None:
        stop_duration = time[-1] - current_stop_start
        if stop_duration >= stop_duration_threshold:
            if current_stop_area == 'safe':
                safe_stops.append(stop_duration)
            else:
                unsafe_stops.append(stop_duration)
                
    return safe_stops, unsafe_stops

def calculate_area_covered_and_speeds(pos_x, pos_z, speeds, time, outer_polygon, inner_polygon):
    """Main function that coordinates all calculations."""
    trajectory_points = [Point(x, z) for x, z in zip(pos_x, pos_z)]

    # Calculate all metrics using separate functions
    coverage = calculate_trajectory_coverage(pos_x, pos_z, outer_polygon, inner_polygon)
    safe_speeds, unsafe_speeds = calculate_speeds_by_area(trajectory_points, speeds, inner_polygon, outer_polygon)
    safe_accelerations, unsafe_accelerations = calculate_accelerations(speeds, time, trajectory_points, inner_polygon, outer_polygon)
    time_in_safe_area, time_in_unsafe_area = calculate_time_in_areas(trajectory_points, time, inner_polygon, outer_polygon)
    safe_stops, unsafe_stops = calculate_stops(trajectory_points, speeds, time, inner_polygon, outer_polygon)

    return (*coverage, safe_speeds, unsafe_speeds, safe_accelerations, unsafe_accelerations,
            time_in_safe_area, time_in_unsafe_area, safe_stops, unsafe_stops)

# functions for t7 and t15 (empty room and dark maze)

def calculate_trajectory_coverage_t7_t15(pos_x, pos_z, outer_polygon):
    """Calculate overall coverage percentage."""
    trajectory_line = LineString(zip(pos_x, pos_z))
    trajectory_area = trajectory_line.buffer(0.1)

    total_intersection = trajectory_area.intersection(outer_polygon)
    total_area = outer_polygon.area
    percentage_covered = (total_intersection.area / total_area) * 100

    return percentage_covered

def calculate_speeds_t7_t15(speeds):
    """Calculate speed metrics."""
    return speeds

def calculate_accelerations_t7_t15(speeds, time):
    """Calculate acceleration metrics."""
    accelerations = []
    
    for i in range(1, len(speeds)):
        
        time_diff = time[i] - time[i-1]
        acceleration = (speeds[i] - speeds[i-1]) / time_diff
        accelerations.append(acceleration)
            
    return accelerations

def calculate_stops_t7_t15(speeds, time, stop_threshold=0.01, stop_duration_threshold=0.5):
    """Calculate all stops in trajectory."""
    stops = []
    current_stop_start = None

    for i, speed in enumerate(speeds):
        if speed < stop_threshold:
            if current_stop_start is None:
                current_stop_start = time[i]
        else:
            if current_stop_start is not None:
                stop_duration = time[i] - current_stop_start
                if stop_duration >= stop_duration_threshold:
                    stops.append(stop_duration)
                current_stop_start = None

    # Check for ongoing stop at the end
    if current_stop_start is not None:
        stop_duration = time[-1] - current_stop_start
        if stop_duration >= stop_duration_threshold:
            stops.append(stop_duration)
                
    return stops

def calculate_trajectory_metrics_t7_t15(pos_x, pos_z, speeds, time, outer_polygon):
    """Main function that coordinates all calculations."""
    # Calculate all metrics using separate functions
    coverage = calculate_trajectory_coverage_t7_t15(pos_x, pos_z, outer_polygon)
    speeds = calculate_speeds_t7_t15(speeds)
    accelerations = calculate_accelerations_t7_t15(speeds, time)
    stops = calculate_stops_t7_t15(speeds, time)

    return coverage, speeds, accelerations, stops