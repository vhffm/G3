"""
Extract Water Mass Fraction, Source Distribution Statistics.
Processing 9000 EJS Outputs Takes ~ 1 Hour on 12 Cores. 
"""

import io_helpers as ioh
import formation_helpers as fh
import multiprocessing as mp
import pandas as pd
import argparse
import sys
import glob
import numpy as np


###############################################################################
# MAIN PROGRAM STARTS HERE
###############################################################################

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--all', action='store_true', \
                    help='Process Full Output Range (9000).')
parser.add_argument('-np', type=int, default=1, \
                    help='Number of Processes')
parser.add_argument('-fout', '--output_file', default='Coordinates_WMF.hdf5', \
                    help='Name of Output File.')
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

# Should Only Be One Directory
if len(dirs) == 1:
    cdir = dirs[0]
else:
    print "!! Script can only process a single output directory."
    print "!! Terminateing."
    sys.exit()

# Define Steps
if args.all:
    nsteps = np.asarray(np.mgrid[0:9e9+1e6:1e6], np.int64)
else:
    nsteps = np.array([ 0, 6e6, 6e7, 6e8, 3e9, 6e9, 9e9], dtype=np.int64)

# The Globbit
globs = glob.glob("%s/Out_*.dat" % cdir)
# Extract run names
# In:  Out_run_03_000057000000
# Out: run_03
run_name = globs[0].strip().split("/")[-1][:-4][4:-13]

# Load Coordinate Outputs
dfo_all = []
for nstep in nsteps:
    fname = "%s/Out_%s_%012d.dat" % (cdir, run_name, nstep)
    dfo = ioh.read_output_and_stack([fname], frame='heliocentric')
    dfo = dfo[dfo.mass < 12.0]
    dfo_all.append(dfo)

# Load Collisions
fname_c = "%s/Collisions_%s.dat" % (cdir, run_name)
dfc = ioh.read_collisions_and_stack([fname_c])

# Compute WMF, Source Regions
targets = []
for istep, nstep in enumerate(nsteps):
    target = [ dfo_all[istep], dfo_all[0], dfc ]
    targets.append(target)

# Loop Steps
print "// Computing WMF & Sources for %i Steps" % len(targets)
pool = mp.Pool(processes=args.np)
result = pool.map(fh.compute_wmf_wrapper, targets, 1)
pool.close()
pool.join()

# Join Data Frame
df = pd.concat(result)
df.reset_index(drop=True, inplace=True)

# Save
print "// Saving to %s" % args.output_file
with pd.HDFStore("%s" % args.output_file, 'w') as store:
    store['df'] = df
