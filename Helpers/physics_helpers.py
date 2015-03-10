"""
Simple Physics Helper Functions.
"""

def rotate_xy(xin, yin, theta):
    """
    Rotate 2d vectors over angle theta.

    @param xin - x-component - [any]
    @param yin - y-component - [any]
    @param theta - rotation angle - [rad]
    @return xout, yout - rotated x/y-components - [any]
    """

    xout = np.cos(theta) * xin - np.sin(theta) * yin
    yout = np.sin(theta) * xin + np.cos(theta) * yin
    return xout, yout


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
