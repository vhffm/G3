"""
Collision Helpers.
"""

import numpy as np
import dynamics_helpers as dh
import vector_helpers as vh
import constants as C


def assign_ell(row):
    """
    Compute Projected Interaction Distance (X-Section) for Collisions.
    Cf. Quintana 2015, Eqn. (4), Fig. 7
    http://adsabs.harvard.edu/abs/2015arXiv151103663Q

    To call this function, use dfc.apply(assign_ell, axis='columns').
    Here, dfc is a DataFrame of Particle Collisions.

    @param: row - Pandas DataFrame Row [Series]
    @return: ell - Projected Interaction Distance (X-Section) [Numpy Float]
    """
    
    B = np.sin(row.theta*C.d2r) * ( row.ri + row.rj )
    
    if (B + row.rj) <= row.ri:
        ell = 2.0 * row.rj
    else:
        ell = row.rj + row.ri - B
    
    return ell


def compute_specific_impact_energy(dfc):
    """
    Compute the Specific Impact Energy for Collisions.
    Cf. Quintana 2015, Eqns. (1,2,3,4,5,6,7)
    http://adsabs.harvard.edu/abs/2015arXiv151103663Q

    The input collision dataframe must have been processed to include
    collisions geometries. We require impact velocities and angle.

    @param: dfc - Genga Collision DataFrame [Pandas DataFrame]
    @return: dfc - Genga Collision DataFrame [Pandas DataFrame]
    """

    # Reduce Mass, Base Specific Impact Energy
    # Eqns. (1,2)
    mu = (dfc.mi * dfc.mj) / (dfc.mi + dfc.mj)           # Msun
    Q = mu * dfc.v_impact**2.0 / 2.0 / (dfc.mi + dfc.mj) # km2/s2
    Q *= 1000.0**2.0                                     # m2/s2
                                                         # = kg m2/s2 1/kg
                                                         # = J/kg

    # Write Uncorrected Q
    dfc['Q'] = Q

    # Compute Geometry Corrected Specific Impact Energy
    # Eqns. (4,5,6,7)
    ell = dfc.apply(assign_ell, axis='columns')
    alpha = ( 3.0 * dfc.rj * ell**2.0 - ell**3.0 ) / ( 4.0 * dfc.rj**3.0 )
    mu_alpha = ( alpha * dfc.mi * dfc.mj ) / ( alpha * dfc.mj + dfc.mi )
    QPrime = Q * mu_alpha / mu

    # Write QPrime
    dfc['QPrime'] = QPrime

    # Compute Shock Wave Corrected Specific Impact Energy
    # Eqn. (3)
    Qs = QPrime * ( 1.0 + dfc.mj/dfc.mi ) * ( 1.0 - np.sin(dfc.theta * C.d2r) )

    # Write to DataFrame
    dfc['Qs'] = Qs

    # Return DataFrame
    return dfc


def reconstruct_geometries(dfc):
    """
    Reconstruct collision geometries by integrating to the boundary.

    Goes through all collisions in a collision file (dataframe), etxracts
    the geometry, integrates to the boundary, and computes impact parameter 
    and impact angle. The impact parameter is normalized to the sum of the
    bodies' radii.

    For collisions dataframe, see io_helpers.py/read_collisions_and_stack().

    We expect certain units in the collision dataframe:
    - Distances are in km.
    - Velocities in km/s.
    - Masses in Earth Masses.

    @param: dfc - Collisions Dataframe [Pandas Dataframe]
    @return: theta_all    - Impact Angles (Degree) [Numpy Float Array]
    @return: b_over_r_all - Normalised Impact Parameter (-) [Numpy Float Array]
    """

    # Compute Relative Position
    dfc['dx'] = dfc.xj - dfc.xi
    dfc['dy'] = dfc.yj - dfc.yi
    dfc['dz'] = dfc.zj - dfc.zi
    dfc['distance'] = np.sqrt(dfc['dx']**2.0 + \
                              dfc['dy']**2.0 + \
                              dfc['dz']**2.0)/(dfc.ri+dfc.rj)

    # Compute Relative Velocity
    dfc['dvx'] = dfc.vxi - dfc.vxj
    dfc['dvy'] = dfc.vyi - dfc.vyj
    dfc['dvz'] = dfc.vzi - dfc.vzj

    # Fix Units
    dfc['dx'] *= C.au2km
    dfc['dy'] *= C.au2km
    dfc['dz'] *= C.au2km
    dfc['ri'] *= C.au2km
    dfc['rj'] *= C.au2km

    theta_all = np.ones(len(dfc)) * np.nan
    b_over_r_all = np.ones(len(dfc)) * np.nan
    v_impact_all = np.ones(len(dfc)) * np.nan   
    for ii, [index, dfc_loc] in enumerate(dfc.iterrows()):

        # Debug. Only do one row.
        # if ii > 3:
        #    break

        # Distance (km)
        x1 = 0.0; y1 = 0.0; z1 = 0.0
        x2 = dfc_loc.dx; y2 = dfc_loc.dy; z2 = dfc_loc.dz

        # Velocity (km/s)
        vx1 = 0.0; vy1 = 0.0; vz1 = 0.0
        vx2 = dfc_loc.dvx; vy2 = dfc_loc.dvy; vz2 = dfc_loc.dvz

        # If vectors are aligned within 180 degree, the dot product
        # is positive. Otherwise, it is negative. For example:
        # vh.dot(2,0,0,1,1,0) ==> 2
        # vh.dot(2,0,0,-1,-1,0) ==> -2
        # We use this to determine whether the impactor is moving towards 
        # or away from the target, i.e. whether we have to integrate forward 
        # or backward.
        if vh.dot(x2, y2, z2, vx2, vy2, vz2) > 0:
            vx2 *= -1.0
            vy2 *= -1.0
            vz2 *= -1.0

        # Mass (kg)
        m1 = dfc_loc.mi * C.mearth
        m2 = dfc_loc.mj * C.mearth

        # Radius (km)
        r1 = dfc_loc.ri
        r2 = dfc_loc.rj

        # Cutoff Distance (km)
        dlo = 0.99*(r1+r2)
        dhi = 1.01*(r1+r2)

        # Estimate Time to Impacts (s), Determine Steps (-)
        dd = np.sqrt(x2**2.0 + y2**2.0 + z2**2.0)
        vv = np.sqrt(vx2**2.0 + vy2**2.0 + vz2**2.0)
        tt = dd/vv
        nsteps = 2048
        dt = tt / float(nsteps)

        # Integrate
        x1_out, y1_out, z1_out, vx1_out, vy1_out, vz1_out, \
        x2_out, y2_out, z2_out, vx2_out, vy2_out, vz2_out, terminated = \
            dh.integrate_system_3d(x1, y1, z1, vx1, vy1, vz1, \
                                   x2, y2, z2, vx2, vy2, vz2, \
                                   m1, m2, dt, 2.5 * nsteps, \
                                   G=C.G, dlo=dlo, dhi=dhi)
            
        # Extract Geometry.
        # Only if the integration managed to put the body on the boundary.
        if terminated:

            # Shift Frame
            vx = np.asarray(vx2_out) - np.asarray(vx1_out)
            vy = np.asarray(vy2_out) - np.asarray(vy1_out)
            vz = np.asarray(vz2_out) - np.asarray(vz1_out)
            x = np.asarray(x2_out) - np.asarray(x1_out)
            y = np.asarray(y2_out) - np.asarray(y1_out)
            z = np.asarray(z2_out) - np.asarray(z1_out)

            # Impact Velocity (km/s)
            v_impact = np.sqrt(vx[-1]**2.0 + vy[-1]**2.0 + vz[-1]**2.0)

            # Extract Geometry
            theta, cos_theta, sin_theta = \
                vh.compute_angle(np.atleast_1d(x[-1]), \
                                 np.atleast_1d(y[-1]), \
                                 np.atleast_1d(z[-1]), \
                                 np.atleast_1d(vx[-1]), \
                                 np.atleast_1d(vy[-1]), \
                                 np.atleast_1d(vz[-1]), \
                                 vanilla=False)
            theta = theta[0]
            b = np.abs(np.sin(theta) * vh.norm(x[-1], y[-1], z[-1]))

            # Shift Angle
            if theta < -np.pi/2.0:
                theta += np.pi
            if theta > np.pi/2.0:
                theta -= np.pi
            theta = np.abs(theta)

            # Add to Master
            theta_all[ii] = theta * C.r2d
            b_over_r_all[ii] = b/(r1+r2)
            v_impact_all[ii] = v_impact

    # Return
    return theta_all, b_over_r_all, v_impact_all
