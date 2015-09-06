"""
Extract Jupiter/Saturn Resonances from a List of Runs. Can be slow.

Computes Secular Resonaces (JS):
- nu_5
- nu_6
- nu_{15}
- nu_{16}

Computes MMR:
- 2:1
- 3:1
- 3:2

With a gas disk present, nu_{16} can appear twice. We store both locations.

Assumes Jupiter/Saturn PIDs 2000/2001.

Usage: python /path/extract_resonances.py -np 2 < dirlist

Here, "dirlist" is a file looking like:
/path/to/run_01
/path/to/run_02
...
/path/to/run_12

Code is embarrassingly parallel. Use the -np to define number of subprocesses.
"""

import sys
import glob
import numpy as np
import resonance_helpers as rh
import io_helpers as ioh
import multiprocessing as mp
import argparse
import constants as C
import pandas as pd


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################

def extract_resonances(cdir):
    """
    Processing Function.
    """

    print "// Processing %s" % cdir

    # Extract run name (again...)
    # In:  Out_run_03_000057000000
    # Out: run_03
    globs = glob.glob("%s/Out_*.dat" % cdir)
    run_name = globs[0].strip().split("/")[-1][:-4][4:-13]
    
    # Running Index for Run
    nrun = int(run_name.strip().split('_')[-1])
    
    # Allocate
    time = np.zeros_like(nsteps) * np.nan
    t_over_tau = np.zeros_like(nsteps) * np.nan
    a_nu_5 = np.zeros_like(nsteps) * np.nan
    a_nu_6 = np.zeros_like(nsteps) * np.nan
    a_nu_15_01 = np.zeros_like(nsteps) * np.nan
    a_nu_15_02 = np.zeros_like(nsteps) * np.nan
    a_nu_16_01 = np.zeros_like(nsteps) * np.nan
    a_nu_16_02 = np.zeros_like(nsteps) * np.nan
    two_to_one = np.zeros_like(nsteps) * np.nan
    three_to_one = np.zeros_like(nsteps) * np.nan
    three_to_two = np.zeros_like(nsteps) * np.nan
    
    # Loop Steps
    for istep, nstep in enumerate(nsteps):
        fname = "%s/Out_%s_%012d.dat" % (cdir, run_name, nstep)
        df = ioh.read_output_and_stack([fname], frame='heliocentric')
        
        # Extract Time & Giant Planets
        time[istep] = df.time.iloc[0]
        t_over_tau[istep] = time[istep]/1.0e6
        a_j = df[df.pid==2000].a.iloc[0]
        m_j = df[df.pid==2000].mass.iloc[0]
        m_j *= C.mearth/C.msun
        a_s = df[df.pid==2001].a.iloc[0]
        m_s = df[df.pid==2001].mass.iloc[0]
        m_s *= C.mearth/C.msun
        
        # Compute Secular Resonance Locations
        a_nu_5_loc, a_nu_6_loc, a_nu_15_loc, a_nu_16_loc = \
            rh.scan_resonances(a_j, a_s, m_j, m_s, t_over_tau[istep])
            
        # Write Back
        a_nu_5[istep] = a_nu_5_loc
        a_nu_6[istep] = a_nu_6_loc
        
        if np.isscalar(a_nu_15_loc):
            a_nu_15_01[istep] = a_nu_15_loc
        else:
            if len(a_nu_15_loc) > 2:
                print "// Found >2 locations for \nu_15. Error?"
            a_nu_15_01[istep] = a_nu_15_loc[0]
            a_nu_15_02[istep] = a_nu_15_loc[1]
            
        if np.isscalar(a_nu_16_loc):
            a_nu_16_01[istep] = a_nu_16_loc
        else:
            if len(a_nu_16_loc) > 2:
                print "// Found >2 locations for \nu_16. Error?"
            a_nu_16_01[istep] = a_nu_16_loc[0]
            a_nu_16_02[istep] = a_nu_16_loc[1]
            
        # Compute Mean Motion Resonance Locations
        two_to_one[istep] = a_j * (2.0/1.0)**(-C.twothirds)
        three_to_one[istep] = a_j * (3.0/1.0)**(-C.twothirds)
        three_to_two[istep] = a_j * (3.0/2.0)**(-C.twothirds)
        
    # Prepare Data Frame
    data = { 'nu_5': a_nu_5, 'nu_6': a_nu_6, \
             'nu_15_01': a_nu_15_01, 'nu_15_02': a_nu_15_02, \
             'nu_16_01': a_nu_16_01, 'nu_16_02': a_nu_16_02, \
             'two_to_one': two_to_one, \
             'three_to_one': three_to_one, \
             'three_to_two': three_to_two, \
             'nrun': np.ones_like(time) * nrun, \
             'time': time, 't_over_tau': t_over_tau, \
             'nstep': np.asarray(nsteps)  }
    
    cols = [ 'nrun', 'time', 't_over_tau', 'nstep', \
             'nu_5', 'nu_6', \
             'nu_15_01', 'nu_15_02', \
             'nu_16_01', 'nu_16_02', \
             'two_to_one', 'three_to_one', 'three_to_two' ]
    
    # Build Data Frame
    df = pd.DataFrame(data, columns=cols)
            
    # Return
    return df


###############################################################################
# MAIN PROGRAM STARTS HERE
###############################################################################

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-np', type=int, default=1, \
                    help='Number of Processes')
args = parser.parse_args()
print "// Using %i Subprocesses" % args.np

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

# Get Steps
nsteps = np.zeros_like(xglobs, dtype=np.int64)
for iglob, xglob in enumerate(sorted(xglobs)):
    # In : /some/dir/Out_run_03_000156000000.dat
    # Out: 156000000
    nsteps[iglob] = int(xglob.strip().split("/")[-1][:-4].split("_")[-1])

# Loop Directories
pool = mp.Pool(processes=args.np)
result = pool.map(extract_resonances, dirs)
pool.close()
pool.join()

# Join Data Frame
df = pd.concat(result)
df.reset_index(drop=True, inplace=True)

# Save
print "// Saving"
with pd.HDFStore('Resonances.hdf5', 'w') as store:
    store['df'] = df

