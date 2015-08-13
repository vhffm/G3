"""
Dynamics Helpers (Integrators).
"""

import numpy as np


def euler_step_3d(x, y, z, vx, vy, vz, ax, ay, az, dt):
    """
    Euler integration in three dimensions. Do one time step.

    *** WARNING ***
    If you don't know about integrators, read up on them. The Euler method 
    (this one) has really bad energy conservation. If you need accuracy or 
    long integration times, use Runge-Kutta or Leapfrog methods. If you're 
    doing Keplerian stuff, consider Wisdom-Holman mappings.

    @param: x, y, z - Positions [Numpy Float (Array)]
    @param: vx, vy, vz - Velocities [Numpy Float (Array)]
    @param: ax, ay, ax - Accelerations [Numpy Float (Array)]
    @return: x, y, z - Positions [Numpy Float (Array)]
    @return: vx, vy, vz - Velocities [Numpy Float (Array)]
    """

    # Velocity
    vx = vx + ax * dt
    vy = vy + ay * dt
    vz = vz + az * dt
    
    # Position
    x = x + vx * dt
    y = y + vy * dt
    z = z + vz * dt

    # Return
    return x, y, z, vx, vy, vz


def integrate_system_3d(x1, y1, z1, vx1, vy1, vz1, \
                        x2, y2, z2, vx2, vy2, vz2, \
                        m1, m2, dt, nsteps, G=1.0, dlo=None, dhi=None):
    """
    Integrate two gravitationally interacting masses.

    The intergration proceeds for "nsteps" steps or until the distance 
    between the objects is
    - (a) in between "dmin" and "dmax" (if both are passed)
    - (b) is larger than "dmin" (if only dmin is passed)
    - (c) is smaller than "dmax" (if only dmax is passed)

    @param: x1, y1, z1    - Body 1 Positions  [Numpy Floats]
    @param: vx1, vy1, vz1 - Body 1 Velocities [Numpy Floats]
    @param: x2, y2, z2    - Body 1 Positions  [Numpy Floats]
    @param: vx2, vy2, vz2 - Body 1 Velocities [Numpy Floats]
    @param: m1, m2        - Masses [Numpy Floats]
    @param: dt            - Time Step [Numpy Float]
    @param: G=1.0         - Gravitational Constant [Numpy Float]
    @param: dhi=None      - Maximum Allowed Separation [Numpy Float]
    @param: dlo=None      - Minimum Allowed Separation [Numpy Float]

    @return: x1_all, y1_all, z1_all    - Body 1 Pos. Time Series [NP Float Arr]
    @return: vx1_all, vy1_all, vz1_all - Body 1 Vel. Time Series [NP Float Arr]
    @return: x2_all, y2_all, z2_all    - Body 2 Pos. Time Series [NP Float Arr]
    @return: vx2_all, vy2_all, vz2_all - Body 2 Vel. Time Series [NP Float Arr]
    @return: terminate                 - Terminated at Boundary? [Boolean]
    """
    
    # Output Arrays - Position
    x1_all = []; y1_all = []; z1_all = []
    x2_all = []; y2_all = []; z2_all = []
    
    # Output Arrays - Velocity
    vx1_all = []; vy1_all = []; vz1_all = []
    vx2_all = []; vy2_all = []; vz2_all = []
    
    # Loop Me
    for nstep in range(int(nsteps)):
        
        # Distance
        dd = (x1-x2)**2.0 + (y1-y2)**2.0 + (z1-z2)**2.0
        
        # if nstep == 0:
        #     print np.sqrt(dd)

        # Forces
        fff = - G * m1 * m2 / dd

        # Accelerations
        ax1 = fff / m1 * (x1-x2)/np.sqrt(dd)
        ax2 = - fff / m2 * (x1-x2)/np.sqrt(dd)

        ay1 = fff / m1 * (y1-y2)/np.sqrt(dd)
        ay2 = - fff / m2 * (y1-y2)/np.sqrt(dd)
        
        az1 = fff / m1 * (z1-z2)/np.sqrt(dd)
        az2 = - fff / m2 * (z1-z2)/np.sqrt(dd)

        # Euler Step
        x1, y1, z1, vx1, vy1, vz1 = \
            euler_step_3d(x1, y1, z1, vx1, vy1, vz1, ax1, ay1, az1, dt)
        x2, y2, z2, vx2, vy2, vz2 = \
            euler_step_3d(x2, y2, z2, vx2, vy2, vz2, ax2, ay2, az2, dt)

        # Append - Position
        x1_all.append(x1)
        x2_all.append(x2)
        y1_all.append(y1)
        y2_all.append(y2)
        z1_all.append(z1)
        z2_all.append(z2)
        
        # Append - Velocity
        vx1_all.append(vx1)
        vx2_all.append(vx2)
        vy1_all.append(vy1)
        vy2_all.append(vy2)
        vz1_all.append(vz1)
        vz2_all.append(vz2)
        
        # Terminate?
        terminate = False
        if (dlo) and (not dhi) and (np.sqrt(dd) > dlo):
            terminate = True
        if (dhi) and (not dlo) and (np.sqrt(dd) < dhi):
            terminate = True
        if (dlo) and (dhi) and (np.sqrt(dd) > dlo) and (np.sqrt(dd) < dhi):
            terminate = True
        if terminate:
            # print "Terminated @ d/dlo, d/dhi = %.4f, %.4f" % \
            #     (np.sqrt(dd)/dlo, np.sqrt(dd)/dhi)
            break
            
    # Throw Notify
    # if not terminate:
    #     print "Failed To Reach Boundary"
        
    # Return
    return x1_all, y1_all, z1_all, vx1_all, vy1_all, vz1_all, \
        x2_all, y2_all, z2_all, vx2_all, vy2_all, vz2_all, \
        terminate
    
