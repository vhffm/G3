"""
Gas Disk Helpers.
"""

import numpy as np
import constants as C


def get_volume_density(r, z=0.0, Sigma0=2000.0):
    """
    Compute Volume Density for Morishima+ 2010 Gas Disk.

    @param: r - Radial Location (AU)
    @param: z - Vertical Location (AU)
    @param: Sigma0 - Surface Density @ 1 AU (g/cm2) (Default 2000)
    @return: rho - Volume Density (kg/m3)
    """

    # {r,z}  = au
    # Sigma0 = g cm-2
    
    # Radius
    r0 = 1.0               # AU
    r0 *= C.au2km * 1000.0 # m
    r *= C.au2km * 1000.0  # m
    
    # Height
    z *= C.au2km * 1000.0 # m
    
    # Gas Stuff
    gamma = 1.4          # -
    T0 = 280             # K
    kB = 1.3806488e-23   # m2 kg s-2 K-1
    mu = 1.660538921e-27 # kg
    R = kB / (2.0 * mu)  # m2 s-2 K-1
    G = 6.67384e-11      # m3 kg-1 s-2
    GM = G * C.msun      # m3 s-2

    # Pulled from Genga (gas.cu)
    F1 = np.sqrt(gamma * R * T0 / GM) # 1/sqrt(m)
    
    # Scale Height
    h = F1 * r0**0.25 * r**1.25     # m
    
    # Surface Density
    Sigma0 /= 1000.0                # kg cm-2
    Sigma0 *= 100.0 * 100.0         # kg m-2
    Sigma = Sigma0 * (r/r0)**(-1.0) # kg m-2
    
    # Density
    rho = (Sigma/np.sqrt(C.twopi)/h) * np.exp(-z**2.0/(2.0 * h**2.0)) # kg m-3
    
    # kg m-3
    return rho


def get_volume_density_mmsn(r, z):
    """
    Compute Volume Density (r,z) for Minimum Mass Solar Nebula (Hayashi 1981).

    @param: r - Radial Location (AU)
    @param: z - Vertical Location (AU)
    @return: rho - Volume Density (kg/m3)
    """

    # AU
    r0 = 1.0

    # g/cm3
    rho = 1.4e-9 * (r/r0)**(-11./4.) * \
          np.exp( - ( ( z / ( 0.047 * (r/r0)**(5./4.) ) )**2.0 ) )

    # kg/m3
    rho /= 1000.0
    rho *= 100.0**3.0

    # Return
    return rho

