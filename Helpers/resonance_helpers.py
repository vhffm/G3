"""
Helpers for Solar System Resonances.

*** Quick Ref -- Secular Resonances, Linear Theory, Two Massive Bodies ***
- g_{1,2}: JS Perigree (omega) precession, eccentricity pumping.
- f_{1,2}: JS Line of nodes (Omega) precession, inclination pumping.
"""

import numpy as np
import constants as C
import potential_helpers as ph


def laplace_coefficient(s, j, alpha):
    """
    Compute Laplace Coefficoents beta_s^j(alpha) from Series Expansion.
    
    Cf. Murray & Dermott (1999), Eq. (6.67), Page 237
        http://adsabs.harvard.edu/abs/1999ssd..book.....M
    Cf. Mardling (2013), App C., Fig. C1
        http://arxiv.org/pdf/1308.0607.pdf
    Cf. Predictability, Stability, and Chaos in N-Body Dynamical Systems,
        Annex A
        https://books.google.ch/books?id=zV_TBwAAQBAJ

    NB: Series convergence is guaranteed for alpha < 1.

    See Murray & Dermott (1999), Eq. (7.33), Page 280 for Test Case.

    @param: s - Half Integer (1/2, 3/2, ...)
    @param: j - Integer (1, 2, 3, ...)
    @param: alpha - Semi-Major Axis Ratio (alpha < 1)
    @return: beta_s^j(alpha) - [Float]
    """
    
    # Leading Term
    if j == 0:
        ss = 1.0
    else:
        ss = s
        for jj in range(1,j):
            ss *= (s+jj)    
    term_01 = ss / np.math.factorial(j) * alpha**j
    
    # Cast
    j = float(j)
    
    # Terms in Square Brackets
    term_02 = s * ( s + j ) / \
                  ( j + 1.0 ) * alpha**2.0 
    
    term_03 = s * ( s + 1.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) \
                  / \
                  ( 2.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) ) * alpha**4.0

    term_04 = s * ( s + 1.0 ) * \
                ( s + 2.0 ) * \
                ( s + j ) * \
                ( s + j + 1.0 ) * \
                ( s + j + 2.0 ) \
                / \
                ( 2.0 * \
                  3.0 * \
                  ( j + 1.0 ) * \
                  ( j + 2.0 ) * \
                  ( j + 3.0 ) ) * alpha**6.0
        
    term_05 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) \
                  / \
                  ( 2.0 * \
                    3.0 * \
                    4.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) ) * alpha**8.0
        
    term_06 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + 4.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) * \
                  ( s + j + 4.0 ) \
                  / \
                  ( \
                    2.0 * \
                    3.0 * \
                    4.0 * \
                    5.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) * \
                    ( j + 5.0 ) ) * alpha**10.0
        
    term_07 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + 4.0 ) * \
                  ( s + 5.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) * \
                  ( s + j + 4.0 ) * \
                  ( s + j + 5.0 ) \
                  / \
                  ( \
                    2.0 * \
                    3.0 * \
                    4.0 * \
                    5.0 * \
                    6.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) * \
                    ( j + 5.0 ) * \
                    ( j + 6.0 ) ) * alpha**12.0
    
    term_08 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + 4.0 ) * \
                  ( s + 5.0 ) * \
                  ( s + 6.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) * \
                  ( s + j + 4.0 ) * \
                  ( s + j + 5.0 ) * \
                  ( s + j + 6.0 ) \
                  / \
                  ( \
                    2.0 * \
                    3.0 * \
                    4.0 * \
                    5.0 * \
                    6.0 * \
                    7.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) * \
                    ( j + 5.0 ) * \
                    ( j + 6.0 ) * \
                    ( j + 7.0 ) ) * alpha**14.0
        
    term_09 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + 4.0 ) * \
                  ( s + 5.0 ) * \
                  ( s + 6.0 ) * \
                  ( s + 7.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) * \
                  ( s + j + 4.0 ) * \
                  ( s + j + 5.0 ) * \
                  ( s + j + 6.0 ) * \
                  ( s + j + 7.0 ) \
                  / \
                  ( \
                    2.0 * \
                    3.0 * \
                    4.0 * \
                    5.0 * \
                    6.0 * \
                    7.0 * \
                    8.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) * \
                    ( j + 5.0 ) * \
                    ( j + 6.0 ) * \
                    ( j + 7.0 ) * \
                    ( j + 8.0 ) ) * alpha**16.0
        
    term_10 = s * ( s + 1.0 ) * \
                  ( s + 2.0 ) * \
                  ( s + 3.0 ) * \
                  ( s + 4.0 ) * \
                  ( s + 5.0 ) * \
                  ( s + 6.0 ) * \
                  ( s + 7.0 ) * \
                  ( s + 8.0 ) * \
                  ( s + j ) * \
                  ( s + j + 1.0 ) * \
                  ( s + j + 2.0 ) * \
                  ( s + j + 3.0 ) * \
                  ( s + j + 4.0 ) * \
                  ( s + j + 5.0 ) * \
                  ( s + j + 6.0 ) * \
                  ( s + j + 7.0 ) * \
                  ( s + j + 8.0 ) \
                  / \
                  ( \
                    2.0 * \
                    3.0 * \
                    4.0 * \
                    5.0 * \
                    6.0 * \
                    7.0 * \
                    8.0 * \
                    9.0 * \
                    ( j + 1.0 ) * \
                    ( j + 2.0 ) * \
                    ( j + 3.0 ) * \
                    ( j + 4.0 ) * \
                    ( j + 5.0 ) * \
                    ( j + 6.0 ) * \
                    ( j + 7.0 ) * \
                    ( j + 8.0 ) * \
                    ( j + 9.0 ) ) * alpha**18.0

    # Laplace Coefficient
    beta = term_01 * ( 1.0 + term_02 + term_03 + term_04 + \
                             term_05 + term_06 + term_07 + \
                             term_08 + term_09 + term_10 )
    beta *= 2.0
    
    # Return
    return beta


def secular_frequencies_planets(a_1, a_2, m_1, m_2):
    """
    Compute secular frequencies g_{1,2} and f_{1,2} for a two planet system.

    g_{1,2}: Perigree (omega) precession frequencies, eccentricity pumping.
    f_{1,2}: Line of nodes (Omega) precession frequencies, inclination pumping.

    *** !!! ***
    This is linear theory. Higher order expansions give (a,e,i) dependence.
    *** !!! ***

    Cf. Murray & Dermott (1999), Sect. 7.3, Pages 270, 280, 281
        http://adsabs.harvard.edu/abs/1999ssd..book.....M
    Cf. Nagasawa (2000), Sect. 2.1, Appendix
        http://adsabs.harvard.edu/abs/2000AJ....119.1480N

    @param: a_1, a_2 - Planet Semi-Major Axis of Planets (AU)
    @param: m_1, m_2 - Planet Masses (Solar Masses)
    @returns: g_1, g_2 - Perigee Precession Frequencies (Rad/Twopi)
    @returns: f_1, f_2 - Line of Nodes Precession Frequencies (Rad/Twopi)
    """

    # Solar System Units
    G = 1.0; M = 1.0

    # Mean Motion
    n_1 = np.sqrt(G * (M + m_1) / a_1**3.0)
    n_2 = np.sqrt(G * (M + m_2) / a_2**3.0)

    # Debug
    # print "n_j = %.6f deg/yr" % (n_j * C.r2d * C.twopi)
    # print "n_s = %.6f deg/yr" % (n_s * C.r2d * C.twopi)

    # What are these called?
    P_12 = a_1 / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 2, a_1 / a_2)
    N_12 = a_1 / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 1, a_1 / a_2)

    # Matrix Elements
    A_11 = 2.0 / n_1 / a_1**2.0 * ( G * m_2 * N_12 )
    A_22 = 2.0 / n_2 / a_2**2.0 * ( G * m_1 * N_12 )
    A_12 = - 2.0 * G * P_12 / a_1 / a_2 * np.sqrt( m_1 * m_2 / n_1 / n_2 )

    B_11 = - 2.0 / n_1 / a_1**2.0 * ( G * m_2 * N_12 )
    B_22 = - 2.0 / n_2 / a_2**2.0 * ( G * m_1 * N_12 )
    B_12 = - 2.0 * G * N_12 / a_1 / a_2 * np.sqrt( m_1 * m_2 / n_1 / n_2 )

    # Frequencies (Rad/Twopi) = Eigenvalues of Matrices A and B
    # Converstion: Multiply (* C.r2d * 3600.0 * C.twopi) [Arcsec/Year]
    g_1 = 0.5 * (A_11 + A_22 - np.sqrt( (A_11 - A_22)**2.0 + 4.0 * A_12**2.0 ))
    g_2 = 0.5 * (A_11 + A_22 + np.sqrt( (A_11 - A_22)**2.0 + 4.0 * A_12**2.0 ))
    f_1 = 0.5 * (B_11 + B_22 + np.sqrt( (B_11 - B_22)**2.0 + 4.0 * B_12**2.0 ))
    f_2 = 0.5 * (B_11 + B_22 - np.sqrt( (B_11 - B_22)**2.0 + 4.0 * B_12**2.0 ))

    # Return Frequencies
    return g_1, g_2, f_1, f_2


def secular_frequencies_planets_gas(a_1, a_2, m_1, m_2, t_over_tau):
    """
    Compute secular frequencies g_{1,2} and f_{1,2} for a two planet system.

    g_{1,2}: Perigree (omega) precession frequencies, eccentricity pumping.
    f_{1,2}: Line of nodes (Omega) precession frequencies, inclination pumping.

    *** !!! ***
    This is linear theory. Higher order expansions give (a,e,i) dependence.
    *** !!! ***

    Cf. Murray & Dermott (1999), Sect. 7.3, Pages 270, 280, 281
        http://adsabs.harvard.edu/abs/1999ssd..book.....M
    Cf. Nagasawa (2000), Sect. 2.1, Appendix
        http://adsabs.harvard.edu/abs/2000AJ....119.1480N

    @param: a_1, a_2 - Planet Semi-Major Axis of Planets (AU)
    @param: m_1, m_2 - Planet Masses (Solar Masses)
    @param: t_over_tau - Exponential Decay Factor for Gas Density (-)
    @returns: g_1, g_2 - Perigee Precession Rate (Rad/Second)
    @returns: f_1, f_2 - Line of Nodes Precession Rate (Rad/Second)
    """

    # Gas Terms (Internal CGS Conversion - Needs AU!)
    S_1, T_1 = ph.s_and_t(a_1, t_over_tau)
    S_2, T_2 = ph.s_and_t(a_2, t_over_tau)

    # CGS 
    G = C.G_cgs; Msun = C.msun * 1000.0
    m_1 *= Msun
    m_2 *= Msun
    a_1 *= C.au2km * 1000.0 * 100.0
    a_2 *= C.au2km * 1000.0 * 100.0

    # What are these called?
    P_12 = a_1 / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 2, a_1 / a_2)
    N_12 = a_1 / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 1, a_1 / a_2)

    # Mean Motion
    n_1 = np.sqrt(G * (Msun + m_1) / a_1**3.0)
    n_2 = np.sqrt(G * (Msun + m_2) / a_2**3.0)

    # Matrix Elements
    A_11 = 2.0 / n_1 / a_1**2.0 * ( G * m_2 * N_12 + T_1 )
    A_22 = 2.0 / n_2 / a_2**2.0 * ( G * m_1 * N_12 + T_2 )
    A_12 = - 2.0 * G * P_12 / a_1 / a_2 * np.sqrt( m_1 * m_2 / n_1 / n_2 )

    B_11 = - 2.0 / n_1 / a_1**2.0 * ( G * m_2 * N_12 + S_1 )
    B_22 = - 2.0 / n_2 / a_2**2.0 * ( G * m_1 * N_12 + S_2 )
    B_12 = - 2.0 * G * N_12 / a_1 / a_2 * np.sqrt( m_1 * m_2 / n_1 / n_2 )

    # Frequencies (Rad/Second) = Eigenvalues of Matrices A and B
    # Conversion: Multiply (C.r2d*3600.0) * (3600.0*365.25*24.0) [Arcsec/Year]
    g_1 = 0.5 * (A_11 + A_22 - np.sqrt( (A_11 - A_22)**2.0 + 4.0 * A_12**2.0 ))
    g_2 = 0.5 * (A_11 + A_22 + np.sqrt( (A_11 - A_22)**2.0 + 4.0 * A_12**2.0 ))
    f_1 = 0.5 * (B_11 + B_22 + np.sqrt( (B_11 - B_22)**2.0 + 4.0 * B_12**2.0 ))
    f_2 = 0.5 * (B_11 + B_22 - np.sqrt( (B_11 - B_22)**2.0 + 4.0 * B_12**2.0 ))

    # Return Frequencies
    return g_1, g_2, f_1, f_2


def secular_frequencies_js_today():
    """
    Compute Present Day Jupiter/Saturn System Frequencies. Test Case.

    @returns: g_1, g_2 - Perigee Precession Frequencies (Rad/Twopi)
    @returns: f_1, f_2 - Line of Nodes Precession Frequencies (Rad/Twopi)
    """

    # Semi-Major Axes
    a_j = 5.202545
    a_s = 9.554841

    # Masses
    m_j = 9.54786e-4
    m_s = 2.85837e-4

    # Compute
    g_1, g_2, f_1, f_2 = secular_frequencies_planets(a_j, a_s, m_j, m_s)

    # Debug
    # For the 9-th order series expansion of the Laplace coefficients, we
    # should get the following frequencies
    # g_1 = 9.637822 10^-4 deg/yr
    # g_2 = 6.100148 10^-3 deg/yr
    # f_1 = 0.000000
    # f_2 = -7.063930 10^-3 deg/yr
    # print "g_1 = %.6f 10^-4 deg/yr" % (g_1 * C.r2d * C.twopi * 1.0e4)
    # print "g_2 = %.6f 10^-3 deg/yr" % (g_2 * C.r2d * C.twopi * 1.0e3)
    # print "f_1 = %.6f" % (f_1 * C.r2d * C.twopi)
    # print "f_2 = %.6f 10^-3 deg/yr" % (f_2 * C.r2d * C.twopi * 1.0e3)

    # Return Frequencies
    return g_1, g_2, f_1, f_2


def secular_frequencies_planets_and_asteroid(a_1, a_2, m_1, m_2, a_ast):
    """
    Compute secular frequencies g and f for an asteroid (massless test
    particles) in the presence of two massive planets.

    @param: a_1, a_2 - Planet Semi-Major Axis of Planets (AU)
    @param: m_1, m_2 - Planet Masses (Solar Masses)
    @param: a_ast - Asteroid Semi-Major Axis (AU)
    @returns: g - Asteroid Perigee Precession Frequency (Rad/Twopi)
    @returns: f - Asteroid Line of Nodes Precession Frequencies (Rad/Twopi)
    """

    # Solar System Units
    G = M = GM = 1.0

    # Mean Motion
    n_ast = np.sqrt(GM/a_ast**3.0)

    # What are these called?
    N_1 = a_ast / 8.0 / a_1**2.0 * laplace_coefficient(1.5, 1, a_ast / a_1)
    N_2 = a_ast / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 1, a_ast / a_2)

    # Frequencies (Rad/Twopi)
    # Converstion: Multiply (* C.r2d * 3600.0 * C.twopi) [Arcsec/Year]
    g = 2.0 * G / n_ast / a_ast**2.0 * ( m_1 * N_1 + m_2 * N_2 )
    # f = - 2.0 * G / n_ast / a_ast**2.0 * ( m_1 * N_1 + m_2 * N_2 )
    f = -g

    # Return Frequencies 
    return g, f


def secular_frequencies_planets_and_asteroid_gas(a_1, a_2, m_1, m_2, \
                                                 a_ast, t_over_tau):
    """
    Compute secular frequencies g and f for an asteroid (massless test
    particles) in the presence of two massive planets.

    Be careful if you request too many (>~10) samples for a_ast. The code
    runs a for-loop (in Fortran) to integrate the disk potential for each
    value of a_ast, which can take a long time if you request many samples.

    @param: a_1, a_2 - Planet Semi-Major Axis of Planets (AU)
    @param: m_1, m_2 - Planet Masses (Solar Masses)
    @param: a_ast - Asteroid Semi-Major Axis (AU)
    @param: t_over_tau - Exponential Decay Factor for Gas Density (-)
    @returns: g - Asteroid Perigee Precession Rate (Rad/Second)
    @returns: f - Asteroid Line of Nodes Precession Frequencies (Rad/Second)
    """

    # Gas Terms (Internal CGS Conversion - Needs AU!)
    if np.isscalar(a_ast):
        S_ast, T_ast = ph.s_and_t(a_ast, t_over_tau)
    else:
        S_ast = np.zeros_like(a_ast) * np.nan
        T_ast = np.zeros_like(a_ast) * np.nan
        for ia, a_ast_loc in enumerate(a_ast):
            S_ast[ia], T_ast[ia] = ph.s_and_t(a_ast_loc, t_over_tau)

    # CGS 
    G = C.G_cgs; Msun = C.msun * 1000.0
    m_1 *= Msun
    m_2 *= Msun
    a_1 *= C.au2km * 1000.0 * 100.0
    a_2 *= C.au2km * 1000.0 * 100.0
    a_ast *= C.au2km * 1000.0 * 100.0

    # Mean Motion
    n_ast = np.sqrt(G*Msun/a_ast**3.0)

    # What are these called?
    N_1 = a_ast / 8.0 / a_1**2.0 * laplace_coefficient(1.5, 1, a_ast / a_1)
    N_2 = a_ast / 8.0 / a_2**2.0 * laplace_coefficient(1.5, 1, a_ast / a_2)

    # Frequencies (Rad/s)
    # Conversion: Multiply (C.r2d*3600.0) * (3600.0*365.25*24.0) [Arcsec/Year]
    g = 2.0 * G / n_ast / a_ast**2.0 * ( m_1 * N_1 + m_2 * N_2 ) + \
        2.0 / n_ast / a_ast**2.0 * T_ast
    f = - 2.0 * G / n_ast / a_ast**2.0 * ( m_1 * N_1 + m_2 * N_2 ) - \
        2.0 / n_ast / a_ast**2.0 * S_ast

    # Return Frequencies 
    return g, f


def secular_resonance_location(a_1, a_2, m_1, m_2):
    """
    Compute Location of Secular Resonances.

    @param: a_1, a_2 - Planet Semi-Major Axis of Planets (AU)
    @param: m_1, m_2 - Planet Masses (Solar Masses)
    @returns: a_nu_5, a_nu_6 - Perigee (omega) Precession Res. (AU)
    @returns: a_nu_15, a_nu_16 - Line of Nodes (Omega) Precession Res. (AU)
    """

    # Scanning Range, Semi-Major Axis
    a_ast = np.mgrid[0.1:5.0:0.01]

    # Compute Planetary Frequencies
    g_1, g_2, f_1, f_2 = secular_frequencies_planets(a_1, a_2, m_1, m_2)

    # Compute Asteroid (Massless Test Particle) Frequencies
    g, f = secular_frequencies_planets_and_asteroid(a_1, a_2, m_1, m_2, a_ast)

    # Locate
    arg_nu_5 = np.argwhere(np.diff(np.sign(g-g_1))**2.0 > 1.0e-16).squeeze()
    arg_nu_6 = np.argwhere(np.diff(np.sign(g-g_2))**2.0 > 1.0e-16).squeeze()
    arg_nu_15 = np.argwhere(np.diff(np.sign(f-f_1))**2.0 > 1.0e-16).squeeze()
    arg_nu_16 = np.argwhere(np.diff(np.sign(f-f_2))**2.0 > 1.0e-16).squeeze()

    a_nu_5  = a_ast[arg_nu_5]
    a_nu_6  = a_ast[arg_nu_6]
    a_nu_15  = a_ast[arg_nu_15]
    a_nu_16  = a_ast[arg_nu_16]

    # Write NaN for Non-Existing Resonances
    if (not np.isscalar(a_nu_5)) and len(a_nu_5) == 0:
        a_nu_5 = np.nan
    if (not np.isscalar(a_nu_6)) and len(a_nu_6) == 0:
        a_nu_6 = np.nan
    if (not np.isscalar(a_nu_15)) and len(a_nu_15) == 0:
        a_nu_15 = np.nan
    if (not np.isscalar(a_nu_16)) and len(a_nu_16) == 0:
        a_nu_16 = np.nan

    # Return
    return a_nu_5, a_nu_6, a_nu_15, a_nu_16


def secular_resonance_location_js_today():
    """
    Compute Present Day JS Secular Resonance Locations. Test Case.

    @returns: a_nu_5, a_nu_6 - Perigee (omega) Precession Res. (AU)
    @returns: a_nu_15, a_nu_16 - Line of Nodes (Omega) Precession Res. (AU)
    """

    # Semi-Major Axes
    a_j = 5.202545
    a_s = 9.554841

    # Masses
    m_j = 9.54786e-4
    m_s = 2.85837e-4

    # Compute
    a_nu_5, a_nu_6, a_nu_15, a_nu_16 = \
        secular_resonance_location(a_j, a_s, m_j, m_s)

    # Debug
    # The resonance locations should be at
    # nu_5  @ a = 0.62 AU
    # nu_6  @ a = 1.84 AU
    # nu_15 @ a = NaN     (DOES NOT EXIST, f_1 = 0)
    # nu_16 @ a = 1.97 AU
    # print "nu_5  @ a = %.2f AU" % a_nu_5
    # print "nu_6  @ a = %.2f AU" % a_nu_6
    # print "nu_15 @ a = %.2f AU" % a_nu_15
    # print "nu_16 @ a = %.2f AU" % a_nu_16

    # Return
    return a_nu_5, a_nu_6, a_nu_15, a_nu_16
