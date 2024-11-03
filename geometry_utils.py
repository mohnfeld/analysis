# geometry_utils.py

import numpy as np
import numpy.linalg as la

def dist(a, b):
    return np.sqrt(np.square(a-b).sum(axis=1))

def minimum_distance(v, w, p):
  l2 = (dist(v, w) ** 2)[0]
  if l2 == 0.0:
    return dist(p, v)
  t = np.clip(np.dot(p - v, np.transpose(w - v, [1,0])) / l2, 0, 1)
  projection = v + t * (w - v)
  return dist(p, projection)

def mean_center_distance(points, center_point):
    dists = dist(points, center_point)
    return np.mean(dists)

def mean_edge_distance(points, corner_points):
    dists = np.zeros((corner_points.shape[0], len(points)))
    rotated_corner_points = np.roll(corner_points, 1)
    point_pairs = [(corner_points[i], rotated_corner_points[i]) for i in range(len(corner_points))]
    for idx, (p1,p2) in enumerate(point_pairs):
        dists[idx,:] = minimum_distance(np.expand_dims(p1,0),
                                        np.expand_dims(p2,0),
                                        points)
    return np.mean(np.min(dists, axis=0))

def mvee(points, tol=0.1):
    """
    Finds the ellipse equation in "center form"
    (x-c).T * A * (x-c) = 1
    """
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T
    err = tol+1.0
    u = np.ones(N)/N
    while err > tol:
        # assert u.sum() == 1 # invariant
        X = np.dot(np.dot(Q, np.diag(u)), Q.T)
        M = np.diag(np.dot(np.dot(Q.T, la.inv(X)), Q))
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
        print(err)
    c = np.dot(u, points)
    A = la.inv(np.dot(np.dot(points.T, np.diag(u)), points)
               - np.multiply.outer(c, c))/d
    return A, c

def get_ellipsoid_axes(points):
    A, centroid = mvee(points)
    U, D, V = la.svd(A)
    # x, y radii.
    rx, ry = 1./np.sqrt(D)
    # Major and minor semi-axis of the ellipse.
    dx, dy = 2 * rx, 2 * ry
    a, b = max(dx, dy), min(dx, dy)
    print('a, b: {}, {}'.format(a, b))
    return a, b

def get_eccentricity(a,b):
    # Eccentricity
    e = np.sqrt(a ** 2 - b ** 2) / a
    return e

def get_path_length(points):
    return np.sum(np.sqrt((np.diff(points, axis=0) ** 2).sum(axis=1)))

def get_focus_and_eccentricity(points):
    a, b = get_ellipsoid_axes(points)
    focus = get_focus(points, a, b)
    eccentricity = get_eccentricity(a, b)
    return focus, eccentricity
def get_focus(points, a, b):
    length = get_path_length(points)
    return 1 - a*b/(length**2 / (4*np.pi))

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

def point_to_line_distance(x, y, x1, y1, x2, y2):
    A = x - x1
    B = y - y1
    C = x2 - x1
    D = y2 - y1

    dot = A * C + B * D
    len_sq = C * C + D * D
    param = -1
    if len_sq != 0:
        param = dot / len_sq

    if param < 0:
        xx = x1
        yy = y1
    elif param > 1:
        xx = x2
        yy = y2
    else:
        xx = x1 + param * C
        yy = y1 + param * D

    dx = x - xx
    dy = y - yy
    return np.sqrt(dx * dx + dy * dy)

def shortest_distance_to_walls(point, coords):
    x, y = point
    min_distance = float('inf')

    for i in range(len(coords)):
        x1, y1 = coords[i]
        x2, y2 = coords[(i + 1) % len(coords)]  # Wrap around to the first point
        distance = point_to_line_distance(x, y, x1, y1, x2, y2)
        min_distance = min(min_distance, distance)

    return min_distance
