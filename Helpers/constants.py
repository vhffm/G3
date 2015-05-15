"""
Constants and Conversions.
"""

import numpy as np

# Masses
msun     = 1.98844e30 # kg
mmercury = 3.28500e23 # kg
mvenus   = 4.86700e24 # kg
mearth   = 5.97219e24 # kg
mmars    = 6.39000e23 # kg
mjupiter = 1.89813e27 # kg
msaturn  = 5.68319e26 # kg
muranus  = 8.68100e25 # kg
mneptune = 1.02400e26 # kg
mpluto   = 1.30900e22 # kg
mmoon    = 7.34767e22 # kg

# Age of the Solar System
t0ss = -4.568 # Gyr

# Cross-Section Areas
Aearth = 6371.0 * np.pi**2.0 # km2

# Angle Conversions
r2d = 180.0 / np.pi
d2r = np.pi / 180.0

# Distance Conversions
au2km = 149597871.0 # km

# True Constants
G = 6.67384e-11      # m3/kg/s2
G = G / 1000.0**3.0  # km3/kg/s2

# Convenience
twopi = 2.0 * np.pi

# Genga Specific Conversions
# Velocity Units [Genga Units <> km/s]
# Test Case: Earth = 1 (Genga) = 30 (km/s)
genga_to_kms = (twopi/365.25) * au2km/24.0/3600.0
kms_to_genga = 1.0/genga_to_kms
