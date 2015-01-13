"""
Various Helper Functions for Vector Math.
"""

import numpy as np

# Roll Our Own Cross-Product
def cross(x1,y1,z1,x2,y2,z2):
    xc = y1 * z2 - z1 * y2
    yc = z1 * x2 - x1 * z2
    zc = x1 * y2 - y1 * x2
    return xc, yc, zc

# Roll Our Own Dot-Product
def dot(x1,y1,z1,x2,y2,z2):
    return x1*x2+y1*y2+z1*z2

# Vector Norm
def norm(x1,y1,z1):
    return np.sqrt(dot(x1,y1,z1,x1,y1,z1))

# Angle Between Vectors
def compute_angle(x1,y1,z1,x2,y2,z2,vanilla=True):
    # Dot Product -> Angle
    dots = dot(x1,y1,z1,x2,y2,z2)
    cos_theta = dots / (norm(x1,y1,z1) * norm(x2,y2,z2))
    theta = np.arccos(cos_theta)
    # Cross Product -> Angle>0 for Z>0; Angle<0 for Z<0
    cx, cy, cz = cross(x1,y1,z1,x2,y2,z2)
    theta[cz<0] = -theta[cz<0]
    # Fix NaN; Kinda Kinky
    # Dot product for almost parallel vectors is
    # -1 (antiparallel) or +1 (parallel)
    nandots = dots[np.isnan(theta)]
    nandots[nandots<0] = np.pi
    nandots[nandots>0] = 0.0
    theta[np.isnan(theta)] = nandots
    # Get Sin (Better For Small Angles, But No Orientation Info)
    sin_theta = norm(cx,cy,cz) / (norm(x1,y1,z1) * norm(x2,y2,z2))
    # Return
    if vanilla:
        return theta
    else:
        return theta, cos_theta, sin_theta
