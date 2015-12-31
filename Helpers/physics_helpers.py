"""
Simple Physics Helper Functions.
"""

import numpy as np


def mass2radius(mass, rho=2.0):
    """
    Convert mass to radius. Assuming sphere.

    @param mass - mass - [kg]
    @param rho - density - [g/cc]
    @return radius - radius - [km]

    Mass in kg; density in g/cc
    Default density is 2 g/cc; earth has 5 g/cc
    """

    rho /= 1000.0 # kg/cc
    rho *= (100.0*1000.0)**3.0 # kg/km3
    V = mass/rho # km3
    r = (V/np.pi * 0.75)**(1./3.) # km
    return r


def radius2mass(radius, rho=2.0):
    """
    Convert radius to mass. Assuming sphere.

    @param radius - radius - [km]
    @param rho - density - [g/cc]
    @return mass - mass - [kg]

    Radius is km; density in g/cc
    Default density is 2 g/cc; earth has 5 g/cc
    """

    rho /= 1000.0 # kg/cc
    rho *= (100.0*1000.0)**3.0 # kg/km3
    mass = 4./3. * np.pi* radius**3.0 * rho # kg
    return mass


def rm2density(radius, mass):
    """
    Compute Density from Radius & Mass.

    @param radius - radius - [km]
    @param mass - mass - [kg]
    @return density - density - [g/cm3]
    """

    volume = 4.0/3.0 * np.pi * radius**3.0 # km3
    density = mass / volume  # kg/km3
    density *= 1000.0        # g/km3
    density /= (1000.0)**3.0 # g/m3
    density /= (100.0)**3.0  # g/cm3
    return density
