"""
Extract Statistics for
- Mass
- Number of Particles

Usage: python /path/extract_stats.py < dirlist

Here, "dirlist" is a file looking like:
/path/to/run_01
/path/to/run_02
...
/path/to/run_12
"""

import sys
import glob
import numpy as np
import pandas as pd
import io_helpers as ioh
import constants as C
import multiprocessing as mp
import argparse
import weighted as wq

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################

def extract_stats(cdir):
    """
    Processing Function.
    """

    # Cutoff
    m_cutoff = 2.0e23 # kg
    m_cutoff /= C.mearth # earth masses

    print "// Processing %s" % cdir

    # Extract run name (again...)
    # In:  Out_run_03_000057000000
    # Out: run_03
    globs = glob.glob("%s/Out_*.dat" % cdir)
    run_name = globs[0].strip().split("/")[-1][:-4][4:-13]
    
    #
    # Allocate
    #
    time = np.ones_like(nsteps) * np.nan

    # Disk Mass
    disk_mass = np.ones_like(time) * np.nan
    disk_mass_above = np.ones_like(time) * np.nan
    disk_mass_below = np.ones_like(time) * np.nan

    # Number of Particles
    npart = np.ones_like(time, dtype=np.int32) * np.nan
    npart_above = np.ones_like(time, dtype=np.int32) * np.nan
    npart_below = np.ones_like(time, dtype=np.int32) * np.nan

    # Particle Mass
    p_mass_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_mass_above_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_above_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_mass_below_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_mass_below_max = np.ones_like(time, dtype=np.int32) * np.nan

    # Semi-Major Axis
    p_a_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_a_above_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_above_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_a_below_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_a_below_max = np.ones_like(time, dtype=np.int32) * np.nan

    # Eccentricity
    p_e_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_e_above_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_above_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_e_below_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_e_below_max = np.ones_like(time, dtype=np.int32) * np.nan

    # Inclination
    p_i_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_i_above_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_above_max = np.ones_like(time, dtype=np.int32) * np.nan

    p_i_below_min = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_q10 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_q25 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_q50 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_q75 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_q90 = np.ones_like(time, dtype=np.int32) * np.nan
    p_i_below_max = np.ones_like(time, dtype=np.int32) * np.nan
    
    # Loop all steps
    for iout, nout in enumerate(nsteps):
        
        # If we cannot load the output, we'll stick to NaN.
        try:
            # Load output into dataframe.
            # Make sure to drop planets <= 12 Earth masses.
            fname = "%s/Out_%s_%012d.dat" % (cdir, run_name, nout)
            df = ioh.read_output(fname, frame="heliocentric")
            df = df[df.mass <= 12.0]

            # Cutoff
            df_above = df[df.mass >= m_cutoff]
            df_below = df[df.mass < m_cutoff]
        
            #
            # Compute Statistics
            #
            time[iout] = df.time.iloc[0]

            # Disk Mass
            disk_mass[iout] = np.sum(df.mass)
            disk_mass_above[iout] = np.sum(df_above.mass)
            disk_mass_below[iout] = np.sum(df_below.mass)

            # Number of Particles
            npart[iout] = len(df)
            npart_above[iout] = len(df_above)
            npart_below[iout] = len(df_below)

            # Particle Mass
            p_mass_min[iout] = df.mass.min
            p_mass_q10[iout] = np.percentile(df.mass, 10)
            p_mass_q25[iout] = np.percentile(df.mass, 25)
            p_mass_q50[iout] = np.percentile(df.mass, 50)
            p_mass_q75[iout] = np.percentile(df.mass, 75)
            p_mass_q90[iout] = np.percentile(df.mass, 90)
            p_mass_max[iout] = df.mass.max

            p_mass_above_min[iout] = df_above.mass.min
            p_mass_above_q10[iout] = np.percentile(df_above.mass, 10)
            p_mass_above_q25[iout] = np.percentile(df_above.mass, 25)
            p_mass_above_q50[iout] = np.percentile(df_above.mass, 50)
            p_mass_above_q75[iout] = np.percentile(df_above.mass, 75)
            p_mass_above_q90[iout] = np.percentile(df_above.mass, 90)
            p_mass_above_max[iout] = df_above.mass.max

            p_mass_below_min[iout] = df_below.mass.min
            p_mass_below_q10[iout] = np.percentile(df_below.mass, 10)
            p_mass_below_q25[iout] = np.percentile(df_below.mass, 25)
            p_mass_below_q50[iout] = np.percentile(df_below.mass, 50)
            p_mass_below_q75[iout] = np.percentile(df_below.mass, 75)
            p_mass_below_q90[iout] = np.percentile(df_below.mass, 90)
            p_mass_below_max[iout] = df_below.mass.max

            # Semi-Major Axis
            p_a_min[iout] = df.a.min
            p_a_q10[iout] = wq.quantile_1D(df.a, df.mass, 0.10)
            p_a_q25[iout] = wq.quantile_1D(df.a, df.mass, 0.25)
            p_a_q50[iout] = wq.quantile_1D(df.a, df.mass, 0.50)
            p_a_q75[iout] = wq.quantile_1D(df.a, df.mass, 0.75)
            p_a_q90[iout] = wq.quantile_1D(df.a, df.mass, 0.90)
            p_a_max[iout] = df.a.max

            p_a_above_min[iout] = df_above.a.min
            p_a_above_q10[iout] = wq.quantile_1D(df_above.a, df_above.mass, 0.10)
            p_a_above_q25[iout] = wq.quantile_1D(df_above.a, df_above.mass, 0.25)
            p_a_above_q50[iout] = wq.quantile_1D(df_above.a, df_above.mass, 0.50)
            p_a_above_q75[iout] = wq.quantile_1D(df_above.a, df_above.mass, 0.75)
            p_a_above_q90[iout] = wq.quantile_1D(df_above.a, df_above.mass, 0.90)
            p_a_above_max[iout] = df_above.a.max

            p_a_below_min[iout] = df_below.a.min
            p_a_below_q10[iout] = wq.quantile_1D(df_below.a, df_below.mass, 0.10)
            p_a_below_q25[iout] = wq.quantile_1D(df_below.a, df_below.mass, 0.25)
            p_a_below_q50[iout] = wq.quantile_1D(df_below.a, df_below.mass, 0.50)
            p_a_below_q75[iout] = wq.quantile_1D(df_below.a, df_below.mass, 0.75)
            p_a_below_q90[iout] = wq.quantile_1D(df_below.a, df_below.mass, 0.90)
            p_a_below_max[iout] = df_below.a.max

            # Eccentricity
            p_e_min[iout] = df.e.min
            p_e_q10[iout] = wq.quantile_1D(df.e, df.mass, 0.10)
            p_e_q25[iout] = wq.quantile_1D(df.e, df.mass, 0.25)
            p_e_q50[iout] = wq.quantile_1D(df.e, df.mass, 0.50)
            p_e_q75[iout] = wq.quantile_1D(df.e, df.mass, 0.75)
            p_e_q90[iout] = wq.quantile_1D(df.e, df.mass, 0.90)
            p_e_max[iout] = df.e.max

            p_e_above_min[iout] = df_above.e.min
            p_e_above_q10[iout] = wq.quantile_1D(df_above.e, df_above.mass, 0.10)
            p_e_above_q25[iout] = wq.quantile_1D(df_above.e, df_above.mass, 0.25)
            p_e_above_q50[iout] = wq.quantile_1D(df_above.e, df_above.mass, 0.50)
            p_e_above_q75[iout] = wq.quantile_1D(df_above.e, df_above.mass, 0.75)
            p_e_above_q90[iout] = wq.quantile_1D(df_above.e, df_above.mass, 0.90)
            p_e_above_max[iout] = df_above.e.max

            p_e_below_min[iout] = df_below.e.min
            p_e_below_q10[iout] = wq.quantile_1D(df_below.e, df_below.mass, 0.10)
            p_e_below_q25[iout] = wq.quantile_1D(df_below.e, df_below.mass, 0.25)
            p_e_below_q50[iout] = wq.quantile_1D(df_below.e, df_below.mass, 0.50)
            p_e_below_q75[iout] = wq.quantile_1D(df_below.e, df_below.mass, 0.75)
            p_e_below_q90[iout] = wq.quantile_1D(df_below.e, df_below.mass, 0.90)
            p_e_below_max[iout] = df_below.e.max

            # Inclination
            p_i_min[iout] = df.i.min
            p_i_q10[iout] = wq.quantile_1D(df.i, df.mass, 0.10)
            p_i_q25[iout] = wq.quantile_1D(df.i, df.mass, 0.25)
            p_i_q50[iout] = wq.quantile_1D(df.i, df.mass, 0.50)
            p_i_q75[iout] = wq.quantile_1D(df.i, df.mass, 0.75)
            p_i_q90[iout] = wq.quantile_1D(df.i, df.mass, 0.90)
            p_i_max[iout] = df.i.max

            p_i_above_min[iout] = df_above.i.min
            p_i_above_q10[iout] = wq.quantile_1D(df_above.i, df_above.mass, 0.10)
            p_i_above_q25[iout] = wq.quantile_1D(df_above.i, df_above.mass, 0.25)
            p_i_above_q50[iout] = wq.quantile_1D(df_above.i, df_above.mass, 0.50)
            p_i_above_q75[iout] = wq.quantile_1D(df_above.i, df_above.mass, 0.75)
            p_i_above_q90[iout] = wq.quantile_1D(df_above.i, df_above.mass, 0.90)
            p_i_above_max[iout] = df_above.i.max

            p_i_below_min[iout] = df_below.i.min
            p_i_below_q10[iout] = wq.quantile_1D(df_below.i, df_below.mass, 0.10)
            p_i_below_q25[iout] = wq.quantile_1D(df_below.i, df_below.mass, 0.25)
            p_i_below_q50[iout] = wq.quantile_1D(df_below.i, df_below.mass, 0.50)
            p_i_below_q75[iout] = wq.quantile_1D(df_below.i, df_below.mass, 0.75)
            p_i_below_q90[iout] = wq.quantile_1D(df_below.i, df_below.mass, 0.90)
            p_i_below_max[iout] = df_below.i.max
            
            # Clean up so we don't run out of memory from too many iterations
            del df
            del df_below
            del df_above
            
        except:
            pass
        
    # Assemble dataframe panel
    data = { 'time': time, \
             'disk_mass': disk_mass, \
             'disk_mass_above': disk_mass_above, \
             'disk_mass_below': disk_mass_below, \
             'npart': npart, \
             'npart_above': npart_above, \
             'npart_below': npart_below, \
             'p_mass_min': p_mass_min, \
             'p_mass_q10': p_mass_q10, \
             'p_mass_q25': p_mass_q25, \
             'p_mass_q50': p_mass_q50, \
             'p_mass_q75': p_mass_q75, \
             'p_mass_q90': p_mass_q90, \
             'p_mass_max': p_mass_max, \
             'p_mass_above_min': p_mass_above_min, \
             'p_mass_above_q10': p_mass_above_q10, \
             'p_mass_above_q25': p_mass_above_q25, \
             'p_mass_above_q50': p_mass_above_q50, \
             'p_mass_above_q75': p_mass_above_q75, \
             'p_mass_above_q90': p_mass_above_q90, \
             'p_mass_above_max': p_mass_above_max, \
             'p_mass_below_min': p_mass_below_min, \
             'p_mass_below_q10': p_mass_below_q10, \
             'p_mass_below_q25': p_mass_below_q25, \
             'p_mass_below_q50': p_mass_below_q50, \
             'p_mass_below_q75': p_mass_below_q75, \
             'p_mass_below_q90': p_mass_below_q90, \
             'p_mass_below_max': p_mass_below_max, \
             'p_a_min': p_a_min, \
             'p_a_q10': p_a_q10, \
             'p_a_q25': p_a_q25, \
             'p_a_q50': p_a_q50, \
             'p_a_q75': p_a_q75, \
             'p_a_q90': p_a_q90, \
             'p_a_max': p_a_max, \
             'p_a_above_min': p_a_above_min, \
             'p_a_above_q10': p_a_above_q10, \
             'p_a_above_q25': p_a_above_q25, \
             'p_a_above_q50': p_a_above_q50, \
             'p_a_above_q75': p_a_above_q75, \
             'p_a_above_q90': p_a_above_q90, \
             'p_a_above_max': p_a_above_max, \
             'p_a_below_min': p_a_below_min, \
             'p_a_below_q10': p_a_below_q10, \
             'p_a_below_q25': p_a_below_q25, \
             'p_a_below_q50': p_a_below_q50, \
             'p_a_below_q75': p_a_below_q75, \
             'p_a_below_q90': p_a_below_q90, \
             'p_a_below_max': p_a_below_max, \
             'p_e_min': p_e_min, \
             'p_e_q10': p_e_q10, \
             'p_e_q25': p_e_q25, \
             'p_e_q50': p_e_q50, \
             'p_e_q75': p_e_q75, \
             'p_e_q90': p_e_q90, \
             'p_e_max': p_e_max, \
             'p_e_above_min': p_e_above_min, \
             'p_e_above_q10': p_e_above_q10, \
             'p_e_above_q25': p_e_above_q25, \
             'p_e_above_q50': p_e_above_q50, \
             'p_e_above_q75': p_e_above_q75, \
             'p_e_above_q90': p_e_above_q90, \
             'p_e_above_max': p_e_above_max, \
             'p_e_below_min': p_e_below_min, \
             'p_e_below_q10': p_e_below_q10, \
             'p_e_below_q25': p_e_below_q25, \
             'p_e_below_q50': p_e_below_q50, \
             'p_e_below_q75': p_e_below_q75, \
             'p_e_below_q90': p_e_below_q90, \
             'p_e_below_max': p_e_below_max, \
             'p_i_min': p_i_min, \
             'p_i_q10': p_i_q10, \
             'p_i_q25': p_i_q25, \
             'p_i_q50': p_i_q50, \
             'p_i_q75': p_i_q75, \
             'p_i_q90': p_i_q90, \
             'p_i_max': p_i_max, \
             'p_i_above_min': p_i_above_min, \
             'p_i_above_q10': p_i_above_q10, \
             'p_i_above_q25': p_i_above_q25, \
             'p_i_above_q50': p_i_above_q50, \
             'p_i_above_q75': p_i_above_q75, \
             'p_i_above_q90': p_i_above_q90, \
             'p_i_above_max': p_i_above_max, \
             'p_i_below_min': p_i_below_min, \
             'p_i_below_q10': p_i_below_q10, \
             'p_i_below_q25': p_i_below_q25, \
             'p_i_below_q50': p_i_below_q50, \
             'p_i_below_q75': p_i_below_q75, \
             'p_i_below_q90': p_i_below_q90, \
             'p_i_below_max': p_i_below_max }

    cols = [ 'time', \
             'disk_mass', 'disk_mass_above', 'disk_mass_below', \
             'npart', 'npart_above', 'npart_below', \
             'p_mass_min', \
             'p_mass_q10', \
             'p_mass_q25', \
             'p_mass_q50', \
             'p_mass_q75', \
             'p_mass_q90', \
             'p_mass_max', \
             'p_mass_above_min', \
             'p_mass_above_q10', \
             'p_mass_above_q25', \
             'p_mass_above_q50', \
             'p_mass_above_q75', \
             'p_mass_above_q90', \
             'p_mass_above_max', \
             'p_mass_below_min', \
             'p_mass_below_q10', \
             'p_mass_below_q25', \
             'p_mass_below_q50', \
             'p_mass_below_q75', \
             'p_mass_below_q90', \
             'p_mass_below_max', \
             'p_a_min', \
             'p_a_q10', \
             'p_a_q25', \
             'p_a_q50', \
             'p_a_q75', \
             'p_a_q90', \
             'p_a_max', \
             'p_a_above_min', \
             'p_a_above_q10', \
             'p_a_above_q25', \
             'p_a_above_q50', \
             'p_a_above_q75', \
             'p_a_above_q90', \
             'p_a_above_max', \
             'p_a_below_min', \
             'p_a_below_q10', \
             'p_a_below_q25', \
             'p_a_below_q50', \
             'p_a_below_q75', \
             'p_a_below_q90', \
             'p_a_below_max', \
             'p_e_min', \
             'p_e_q10', \
             'p_e_q25', \
             'p_e_q50', \
             'p_e_q75', \
             'p_e_q90', \
             'p_e_max', \
             'p_e_above_min', \
             'p_e_above_q10', \
             'p_e_above_q25', \
             'p_e_above_q50', \
             'p_e_above_q75', \
             'p_e_above_q90', \
             'p_e_above_max', \
             'p_e_below_min', \
             'p_e_below_q10', \
             'p_e_below_q25', \
             'p_e_below_q50', \
             'p_e_below_q75', \
             'p_e_below_q90', \
             'p_e_below_max', \
             'p_i_min', \
             'p_i_q10', \
             'p_i_q25', \
             'p_i_q50', \
             'p_i_q75', \
             'p_i_q90', \
             'p_i_max', \
             'p_i_above_min', \
             'p_i_above_q10', \
             'p_i_above_q25', \
             'p_i_above_q50', \
             'p_i_above_q75', \
             'p_i_above_q90', \
             'p_i_above_max', \
             'p_i_below_min', \
             'p_i_below_q10', \
             'p_i_below_q25', \
             'p_i_below_q50', \
             'p_i_below_q75', \
             'p_i_below_q90', \
             'p_i_below_max' ]

    # Return
    return pd.DataFrame(data, columns = cols)

###############################################################################
# MAIN PROGRAM STARTS HERE
###############################################################################

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-np', type=int, default=1, \
                    help='Number of Processes')
args = parser.parse_args()

# List of Directories
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    for line in lines:
        dirs.append(line)
    print "// Reading %i Directories" % len(dirs)

# The Globbit
run_names = []
for idir, cdir in enumerate(dirs):
    globs = glob.glob("%s/Out_*.dat" % cdir)
    # Extract run names
    # In:  Out_run_03_000057000000
    # Out: run_03
    run_names.append(globs[0].strip().split("/")[-1][:-4][4:-13])
    # Extract directory w/ most outputs
    if idir == 0:
        lmax = len(globs)
        xglobs = globs
    else:
        # Actual
        if len(globs) > lmax:
            lmax = len(globs)
            xglobs = globs

# Get steps
nsteps = np.zeros_like(xglobs, dtype=np.int64)
for iglob, xglob in enumerate(sorted(xglobs)):
    # In : /some/dir/Out_run_03_000156000000.dat
    # Out: 156000000
    nsteps[iglob] = int(xglob.strip().split("/")[-1][:-4].split("_")[-1])

# Loop directories. Serial or parallel.
df_sts_runs = {}
if args.np == 1:
    for idir, cdir in enumerate(dirs):
        df_sts_runs["%s" % run_names[idir]] = extract_stats(cdir)
else:
    pool = mp.Pool(processes=args.np)
    result = pool.map(extract_stats, dirs)
    pool.close()
    pool.join()
    for idir, cdir in enumerate(dirs):
        df_sts_runs["%s" % run_names[idir]] = result[idir]

# Panel
wp = pd.Panel(df_sts_runs)

# Compute Various Quantiles
print "// Computing Quantiles"
dict_min = {}
dict_q10 = {}
dict_q25 = {}
dict_q50 = {}
dict_q75 = {}
dict_q90 = {}
dict_max = {}
# @todo: Parallelize?
for icolumn, column in enumerate(wp[wp.items[0]].columns.values):
    print "   Column %03d/%03d (%s)" % \
        (icolumn+1, len(wp[wp.items[0]].columns.values), column)
    # wp[wp.items[0]] == wp["run_01"]
    # NB: wp (pandas panel) is organized as 
    #     - wp.items (run_01...run_xx)
    #     - wp.major_axis (time...npart)
    #     - wp.minor_axis (1...XX)
    dict_min[column] = wp.xs(column, axis=2).quantile(q=0.00, axis=1)
    dict_q10[column] = wp.xs(column, axis=2).quantile(q=0.10, axis=1)
    dict_q25[column] = wp.xs(column, axis=2).quantile(q=0.25, axis=1)
    # dict_q50 == median
    dict_q50[column] = wp.xs(column, axis=2).quantile(q=0.50, axis=1)
    dict_q75[column] = wp.xs(column, axis=2).quantile(q=0.75, axis=1)
    dict_q90[column] = wp.xs(column, axis=2).quantile(q=0.90, axis=1)
    dict_max[column] = wp.xs(column, axis=2).quantile(q=1.00, axis=1)
    
df_min = pd.DataFrame(dict_min, columns=wp[wp.items[0]].columns.values)
df_q10 = pd.DataFrame(dict_q10, columns=wp[wp.items[0]].columns.values)
df_q25 = pd.DataFrame(dict_q25, columns=wp[wp.items[0]].columns.values)
# df_q50 == wp.median()
df_q50 = pd.DataFrame(dict_q50, columns=wp[wp.items[0]].columns.values)
df_q75 = pd.DataFrame(dict_q75, columns=wp[wp.items[0]].columns.values)
df_q90 = pd.DataFrame(dict_q90, columns=wp[wp.items[0]].columns.values)
df_max = pd.DataFrame(dict_max, columns=wp[wp.items[0]].columns.values)

# Save
print "// Saving"
with pd.HDFStore("Stats.hdf5", "w") as store:
    store['wp'] = wp
    store['df_min'] = df_min
    store['df_q10'] = df_q10
    store['df_q25'] = df_q25
    store['df_q50'] = df_q50
    store['df_q75'] = df_q75
    store['df_q90'] = df_q90
    store['df_max'] = df_max
