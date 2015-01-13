"""
Various Helper Function for Kepler Orbits.

They usually break down for inc=0 and/or ecc=0.
Might need to upgrade to Regular Equations some day
cf. http://server.faia.upm.es/moda/curso1112/kepler.pdf
"""

import numpy as np
import vector_helpers as vh

def helio2bary(r_h, v_h, m, return_sun=False, central_mass=1.0):
    """
    Convert Helio- to Barycentric Coordinates.
    """

    total_mass = np.sum(m) + central_mass

    # Barycenter Offsets
    r_b = np.sum(r_h * m) / total_mass
    v_b = np.sum(v_h * m) / total_mass

    # New Coordinates
    r = r_h - r_b
    v = v_h - v_b

    if return_sun:
        return r, v, -r_b, -v_b
    else:
        return r, v

def cart2kepX(x, y, z, vx, vy, vz, mass, central_mass=1.0):
    """
    Vectorized version of cart2kep. Much love for many particles.

    @params
    r - (x,y,z) Cartesian Positions
    v - (vx,vy,vz) Cartesian Velocities
    mass - Particle Mass
    central_mass - Mass of Central Object

    @returns
    a - Semi-Major Axis
    ecc - Eccentricity
    inc - Inclination
    Omega - Longitude of the Ascending Node
    omega - Argument of Periapsis
    M - Mean Anomaly at Epoch

    Cf. http://www.bruce-shapiro.com/pair/ElementConversionRecipes.pdf
    """

    # Gravitational Parameter
    G = 1.0
    mu = G * ( central_mass + mass )

    # Angular Momentum Vector
    hx, hy, hz = vh.cross(x, y, z, vx, vy, vz)

    # Laplace-Runge-Lenz Vector
    # Scalar Eccentricity
    tmp_x, tmp_y, tmp_z = vh.cross(vx, vy, vz, hx, hy, hz)
    r_norm = vh.norm(x, y, z)
    lrl_x = tmp_x / mu - x / r_norm
    lrl_y = tmp_y / mu - y / r_norm
    lrl_z = tmp_z / mu - z / r_norm
    ecc = vh.norm(lrl_x, lrl_y, lrl_z)

    # Semi-Major Axis
    a = vh.dot(hx, hy, hz, hx, hy, hz) / ( mu * ( 1.0 - ecc**2.0 ) )

    # Inclination
    h_norm = vh.norm(hx, hy, hz)
    inc = np.arccos(vh.dot(0.0, 0.0, 1.0, hx, hy, hz) / h_norm)

    # Longitude of the Ascending Node
    nx, ny, nz = vh.cross(0.0, 0.0, 1.0, hx, hy, hz)
    n_norm = vh.norm(nx, ny, nz)
    tmp0 = 0.0
    tmpX = np.arccos(vh.dot(1.0, 0.0, 0.0, nx, ny, nz) / n_norm)
    Omega = np.where(inc==0.0, tmp0, tmpX)
    Omega[vh.dot(nx, ny, nz, 0.0, 1.0, 0.0) < 0.0] = \
        2.0 * np.pi - Omega[vh.dot(nx, ny, nz, 0.0, 1.0, 0.0) < 0.0]

    # Argument of Perigee
    # For Zero Inclination, Fall Back to 2D Case
    # http://en.wikipedia.org/wiki/Argument_of_periapsis
    tmp0 = np.arctan2(lrl_y/ecc, lrl_x/ecc)
    tmpX = np.arccos(vh.dot(nx, ny, nz, lrl_x, lrl_y, lrl_z) / (n_norm * ecc))
    omega = np.where(inc==0.0, tmp0, tmpX)
    omega[vh.dot(lrl_x, lrl_y, lrl_z, 0.0, 0.0, 1.0) < 0.0] = \
        2.0 * np.pi - omega[vh.dot(lrl_x, lrl_y, lrl_z, 0.0, 0.0, 1.0) < 0.0]

    # True Anomaly
    theta = np.arccos(vh.dot(lrl_x, lrl_y, lrl_z, x, y, z) / (ecc * r_norm))
    theta[vh.dot(x, y, z, vx, vy, vz) < 0.0] = \
        2.0 * np.pi - theta[vh.dot(x, y, z, vx, vy, vz) < 0.0]

    # Eccentric Anomaly
    E = np.arccos((ecc + np.cos(theta)) / (1 + ecc * np.cos(theta)))
    E[np.logical_and(np.pi < theta, theta < 2.0 * np.pi)] = \
        2.0 * np.pi - E[np.logical_and(np.pi < theta, theta < 2.0 * np.pi)]

    # Mean Anomaly
    M = E - ecc * np.sin(E)

    # Return Set
    return a, ecc, inc, Omega, omega, M

def cart2kep(r, v, mass, central_mass=1.0):
    """
    @params
    r - (x,y,z) Cartesian Positions
    v - (vx,vy,vz) Cartesian Velocities
    mass - Particle Mass
    central_mass - Mass of Central Object

    @returns
    a - Semi-Major Axis
    ecc - Eccentricity
    inc - Inclination
    Omega - Longitude of the Ascending Node
    omega - Argument of Periapsis
    M - Mean Anomaly at Epoch

    Cf. http://www.bruce-shapiro.com/pair/ElementConversionRecipes.pdf
    """

    # Gravitational Parameter
    G = 1.0
    mu = G * ( central_mass + mass )

    # Cartesian Unit Vectors
    iHat = np.array([1., 0., 0.])
    jHat = np.array([0., 1., 0.])
    kHat = np.array([0., 0., 1.])

    # Eccentricity
    h = np.cross(r, v)
    evec = 1. / mu * np.cross(v, h) - r / np.linalg.norm(r)
    ecc = np.linalg.norm(evec)

    # Semi Major Axis
    a = np.dot(h,h) / ( mu * ( 1. - ecc**2. ))

    # Inclination
    inc = np.arccos(np.dot(kHat, h) / np.linalg.norm(h))

    # Longitude of the Ascending Node
    n = np.cross(kHat, h)
    if inc == 0.0:
        Omega = 0.0
    else:
        Omega = np.arccos(np.dot(iHat, n) / np.linalg.norm(n))
        if np.dot(n, jHat) < 0:
            Omega = 2. * np.pi - Omega

    # Argument of Perigee
    # For Zero Inclination, Fall Back to 2D Case
    # http://en.wikipedia.org/wiki/Argument_of_periapsis
    if inc == 0.0:
        omega = np.arctan2(evec[1]/ecc,evec[0]/ecc)
    else:
        omega = np.arccos(np.dot(n, evec) / (np.linalg.norm(n) * ecc))
        if np.dot(evec, kHat) < 0:
            omega = 2. * np.pi - omega

    # True Anomaly
    theta = np.arccos(np.dot(evec, r) / (ecc * np.linalg.norm(r)))
    if np.dot(r,v) < 0:
        theta = 2. * np.pi - theta

    # Eccentric Anomaly
    E = np.arccos((ecc + np.cos(theta)) / (1 + ecc * np.cos(theta)))
    if np.pi < theta and theta < 2. * np.pi:
        E = 2. * np.pi - E

    # Mean Anomaly
    M = E - ecc * np.sin(E)

    # Return Set
    return a, ecc, inc, Omega, omega, M

def kep2cart(a, ecc, inc, Omega, omega, M, mass, central_mass=1.0):
    """
    @params
    a - Semi-Major Axis
    ecc - Eccentricity
    inc - Inclination
    Omega - Longitude of the Ascending Node
    omega - Argument of Periapsis
    M - Mean Anomaly at Epoch

    @returns
    r - (x,y,z) Cartesian Positions
    v - (vx,vy,vz) Cartesian Velocities

    Cf. http://www.bruce-shapiro.com/pair/ElementConversionRecipes.pdf
    """

    # Zero Inclination has no Argument of Perigee
    # Set =0 by Convention
    if inc == 0.0:
        Omega = 0.0

    # Mean Anomaly -> Eccentric Anomaly
    E = nr(M, ecc)

    # PQW Unit Vectors
    Px = np.cos(omega) * np.cos(Omega) - \
         np.sin(omega) * np.cos(inc) * np.sin(Omega)
    Py = np.cos(omega) * np.sin(Omega) + \
         np.sin(omega) * np.cos(inc) * np.cos(Omega)
    Pz = np.sin(omega) * np.sin(inc)

    Qx = - np.sin(omega) * np.cos(Omega) - \
           np.cos(omega) * np.cos(inc) * np.sin(Omega)
    Qy = - np.sin(omega) * np.sin(Omega) + \
           np.cos(omega) * np.cos(inc) * np.cos(Omega)
    Qz =   np.sin(inc) * np.cos(omega)

    P = np.array([Px, Py, Pz])
    Q = np.array([Qx, Qy, Qz])

    # Position
    x = a * (np.cos(E) - ecc) * P + \
        a * np.sqrt(1 - ecc**2.) * np.sin(E) * Q

    # Velocity
    G = 1.0
    mu = G * ( central_mass + mass )
    Edot = np.sqrt(mu / a**3.) / ( 1. - ecc * np.cos(E) )
    v = - a * np.sin(E) * Edot * P + \
          a * np.sqrt(1. - ecc**2.) * np.cos(E) * Edot * Q

    # Return
    return x, v

def nr(M, ecc, epsilon_target=1.0e-5):
    """
    Newton-Raphson Iteration to Compute Eccentric Anomaly from Mean Anomaly.
    """

    Ei = M; ii = 1
    # print "Running Newton-Raphson for M=%.2e, e=%.2f" % ( M, ecc )
    while True:
        Ei1 = Ei - ( Ei - ecc * np.sin(Ei) - M ) / (  1 - ecc * np.cos(Ei)  )
        epsilon = np.abs(Ei1 - Ei)
        Ei = Ei1
        # print "Iteration %i, Residual %.2e" % ( ii, epsilon )
        if epsilon < epsilon_target:
            break
        ii += 1
    # print "Found E=%.2f" % Ei1
    return Ei1

def compute_ellipseX(a, ecc, inc, Omega, omega):
    """
    Compute XYZ Sequence for a Kepler Ellipse.
    Vectorized Version. Expects 1D Arrays Passed.
    """

    # Some Reshaping
    a = a[:,np.newaxis]
    ecc = ecc[:,np.newaxis]
    inc = inc[:,np.newaxis]
    Omega = Omega[:,np.newaxis]
    omega = omega[:,np.newaxis]

    # Eccentric Anomaly
    E = np.linspace(0.0, 2.*np.pi, 128)
    E = E[np.newaxis,:]

    # PQW Unit Vectors
    Px = np.cos(omega) * np.cos(Omega) - \
         np.sin(omega) * np.cos(inc) * np.sin(Omega)
    Py = np.cos(omega) * np.sin(Omega) + \
         np.sin(omega) * np.cos(inc) * np.cos(Omega)
    Pz = np.sin(omega) * np.sin(inc)

    Qx = - np.sin(omega) * np.cos(Omega) - \
           np.cos(omega) * np.cos(inc) * np.sin(Omega)
    Qy = - np.sin(omega) * np.sin(Omega) + \
           np.cos(omega) * np.cos(inc) * np.cos(Omega)
    Qz =   np.sin(inc) * np.cos(omega)

    # Compute Ellipse
    x = a * (np.cos(E) - ecc) * Px + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qx
    y = a * (np.cos(E) - ecc) * Py + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qy
    z = a * (np.cos(E) - ecc) * Pz + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qz

    # Return
    return x, y, z

def compute_ellipse(a, ecc, inc, Omega, omega):
    """
    Compute XYZ Sequence for a Kepler Ellipse.
    """

    # Eccentric Anomaly
    E = np.linspace(0.0, 2.*np.pi, 128)

    # PQW Unit Vectors
    Px = np.cos(omega) * np.cos(Omega) - \
         np.sin(omega) * np.cos(inc) * np.sin(Omega)
    Py = np.cos(omega) * np.sin(Omega) + \
         np.sin(omega) * np.cos(inc) * np.cos(Omega)
    Pz = np.sin(omega) * np.sin(inc)

    Qx = - np.sin(omega) * np.cos(Omega) - \
           np.cos(omega) * np.cos(inc) * np.sin(Omega)
    Qy = - np.sin(omega) * np.sin(Omega) + \
           np.cos(omega) * np.cos(inc) * np.cos(Omega)
    Qz =   np.sin(inc) * np.cos(omega)

    # Compute Ellipse
    x = a * (np.cos(E) - ecc) * Px + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qx
    y = a * (np.cos(E) - ecc) * Py + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qy
    z = a * (np.cos(E) - ecc) * Pz + \
        a * np.sqrt(1.0 - ecc**2.) * np.sin(E) * Qz

    # Return
    return x, y, z

def kep2del(a, e, i, Omega, omega, M):
    """
    Keper to Delaunay Elements.
    Cf. (a) Joachim's Unpublished Paper
        (b) http://www.bourbaphy.fr/chenciner.pdf
    """

    Msolar = 1.99e30 # kg
    m = 5.0 * 5.97e24 / 2048.0
    mu = m**2.0 * Msolar
    Lambda = np.sqrt(mu * a)
    L = np.sqrt(mu * a * ( 1.0 - e**2.0 ))
    Lz = L * np.cos(i)
    return Lambda, L, Lz, Omega, omega, M

def del2poi(Lambda, L, Lz, Omega, omega, M):
    """
    Delaunay to Poincare Elements.
    Cf. (a) Joachim's Unpublished Paper
        (b) http://www.bourbaphy.fr/chenciner.pdf
    """

    lambda_small = M + omega + Omega
    lambda_small = lambda_small % ( 2.0 * np.pi )
    xi_real = np.sqrt(2.0 * ( Lambda - L )) * np.cos(omega + Omega)
    xi_imag = np.sqrt(2.0 * ( Lambda - L )) * np.sin(omega + Omega)
    eta_real = np.sqrt(2.0 * ( L - Lz )) * np.cos(Omega)
    eta_imag = np.sqrt(2.0 * ( L - Lz )) * np.sin(Omega)
    return Lambda, lambda_small, xi_real, xi_imag, eta_real, eta_imag

def kep2metric(a1, a2, e1, e2, i1, i2, \
               Omega1, Omega2, omega1, omega2):
    """
    Computes Metric Distance from Keplerian Elements.
    Cf. http://adsabs.harvard.edu/abs/2008CeMDA.100..169K
    Eq. 21

    This breaks down for small inclinations, eccentricities.
    Here, small changes in coordinates cause omega, Omega to vary hugely.
    Thus, we cannot compute small differences anymore.
     
    Calculate directly from state vector (x,v) in this case.
    See below in cart2metic, cart2metricX.
    """

    # Some sort of normalization
    L = 1.0
    # Help Me 01
    cos_xi = np.cos(i1) * np.cos(i2) + np.sin(i1) * np.sin(i2) * np.cos(Omega1 - Omega2)
    # Help Me 02
    cos_eta = ( np.cos(omega1) * np.cos(omega2) + np.cos(i1) * np.cos(i2) * np.sin(omega1) * np.sin(omega2) ) * np.cos(Omega1 - Omega2) \
            + ( np.cos(i2) * np.cos(omega1) * np.sin(omega2) - np.cos(i1) * np.sin(omega1) * np.cos(omega2) ) * np.sin(Omega1 - Omega2) \
            + np.sin(i1) * np.sin(i2) * np.sin(omega1) * np.sin(omega2)
    # One More Conversion
    p1 = a1 * ( 1.0 - e1**2.0 )
    p2 = a2 * ( 1.0 - e2**2.0 )
    # Hit Me
    rho2 = 1.0 / L * ( p1 + p2 - 2.0 * np.sqrt(p1*p2) * cos_xi ) + ( e1**2.0 + e2**2.0 - 2.0 * e1 * e2 * cos_eta )
    return rho2

def cart2metric(r1, v1, r2, v2):
    """
    Computes Metric Distance.
    Cf. Kholshevnikov 2007, Sec. 3, Eq. (4)
    """
    # mu = G * ( central_mass + mass )
    mu = 1.0
    L = 1.0
    # Specific Angular Momentum Vectors
    h1 = np.cross(r1,v1)
    h2 = np.cross(r2,v2)
    # Laplace-Runge-Lenz Vectors / Eccentricity Vectors
    e1 = np.cross(v1,h1)/mu - r1/np.linalg.norm(r1)
    e2 = np.cross(v2,h2)/mu - r2/np.linalg.norm(r2)
    # Distance (Squared)
    dh = h1 - h2
    de = e1 - e2
    rho2 = 1.0/mu/L * np.dot(dh,dh) + np.dot(de,de)
    return rho2

def cart2metricX(x1, y1, z1, vx1, vy1, vz1, \
                 x2, y2, z2, vx2, vy2, vz2, \
                 vanilla=True):
    """
    Computes Metric Distance.
    Cf. Kholshevnikov 2007, Sec. 3, Eq. (4)
    """
    # Gravitational Constant
    mu = 1.0
    # Scale Factors
    L = 1.0
    L1 = 1.0
    # Specific Angular Momentum Vectors
    hx1, hy1, hz1 = vh.cross(x1, y1, z1, vx1, vy1, vz1)
    hx2, hy2, hz2 = vh.cross(x2, y2, z2, vx2, vy2, vz2)
    # Laplace-Runge-Lenz Vectors / Eccentricity Vectors
    vx1, vy1, vz1 = vh.cross(vx1, vy1, vz1, hx1, hy1, hz1)
    vx2, vy2, vz2 = vh.cross(vx2, vy2, vz2, hx2, hy2, hz2)
    ex1 = vx1 / mu - x1 / vh.norm(x1, y1, z1)
    ey1 = vy1 / mu - y1 / vh.norm(x1, y1, z1)
    ez1 = vz1 / mu - z1 / vh.norm(x1, y1, z1)
    ex2 = vx2 / mu - x2 / vh.norm(x2, y2, z2)
    ey2 = vy2 / mu - y2 / vh.norm(x2, y2, z2)
    ez2 = vz2 / mu - z2 / vh.norm(x2, y2, z2)
    # Energy Constants
    E1 = vh.dot(vx1, vy1, vz1, vx1, vy1, vz1) / 2.0 - mu / vh.norm(x1, y1, z1)
    E2 = vh.dot(vx2, vy2, vz2, vx2, vy2, vz2) / 2.0 - mu / vh.norm(x2, y2, z2)
    # Coordinate Distances (Squared)
    dh2 = (hx1 - hx2)**2.0 + (hy1 - hy2)**2.0 + (hz1 - hz2)**2.0
    de2 = (ex1 - ex2)**2.0 + (ey1 - ey2)**2.0 + (ez1 - ez2)**2.0
    dE2 = (E1 - E2)**2.0
    # Total Distances (Squared)
    rho2 = 1.0/mu/L * dh2 + de2
    rho2e = 1.0/mu/L * dh2 + de2 + L1**2.0 / mu**2.0 * dE2
    # Return Values
    if vanilla:
        return rho2
    else:
        return rho2, rho2e, dh2, de2, dE2
