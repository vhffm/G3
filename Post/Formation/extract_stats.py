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

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################

def extract_stats(cdir):
    """
    Processing Function.
    """

    print "// Processing %s" % cdir

    # Extract run name (again...)
    # In:  Out_run_03_000057000000
    # Out: run_03
    globs = glob.glob("%s/Out_*.dat" % cdir)
    run_name = globs[0].strip().split("/")[-1][:-4][4:-13]
    
    # Allocate
    time = np.ones_like(nsteps) * np.nan
    mass_total = np.ones_like(time) * np.nan
    mass_above_cutoff = np.ones_like(time) * np.nan
    mass_below_cutoff = np.ones_like(time) * np.nan
    npart_total = np.ones_like(time, dtype=np.int32) * np.nan
    npart_above_cutoff = np.ones_like(time, dtype=np.int32) * np.nan
    npart_below_cutoff = np.ones_like(time, dtype=np.int32) * np.nan
    
    # Loop all steps
    for iout, nout in enumerate(nsteps):
        
        # If we cannot load the output, we'll stick to NaN.
        try:
            # Load output into dataframe.
            # Make sure to drop planets >12 Earth masses.
            fname = "%s/Out_%s_%012d.dat" % (cdir, run_name, nout)
            df_out_loc = ioh.read_output(fname, frame="heliocentric")
            df_out_loc = df_out_loc[df_out_loc.mass <= 12.0]
        
            # Compute statistics
            time[iout] = df_out_loc.time.iloc[0]
            mass_total[iout] = np.sum(df_out_loc.mass)
            mass_above_cutoff[iout] = \
                np.sum(df_out_loc[df_out_loc.mass >= m_cutoff].mass)
            mass_below_cutoff[iout] = \
                np.sum(df_out_loc[df_out_loc.mass < m_cutoff].mass)
            npart_total[iout] = len(df_out_loc)
            npart_above_cutoff[iout] = np.sum(df_out_loc.mass >= m_cutoff)
            npart_below_cutoff[iout] = np.sum(df_out_loc.mass < m_cutoff)
            
            # Clean up so we don't run out of memory from too many iterations
            del df_out_loc
            
        except:
            pass
        
    # Assemble dataframe panel
    data = { 'time': time, \
             'mass_total': mass_total, \
             'mass_above_cutoff': mass_above_cutoff, \
             'mass_below_cutoff': mass_below_cutoff, \
             'npart_total': npart_total, \
             'npart_above_cutoff': npart_above_cutoff, \
             'npart_below_cutoff': npart_below_cutoff }
    cols = [ 'time', \
             'mass_total', 'mass_above_cutoff', 'mass_below_cutoff', \
             'npart_total', 'npart_above_cutoff', 'npart_below_cutoff']

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

# Cutoff
m_cutoff = 2.0e23 # kg
m_cutoff /= C.mearth # earth masses

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

# Compute Medians over Runs
print "// Computing Medians"
df_medians = wp.median(axis=0)

# Save
print "// Saving"
with pd.HDFStore("Stats.hdf5", "w") as store:
    store['wp'] = wp
    store['df_medians'] = df_medians
