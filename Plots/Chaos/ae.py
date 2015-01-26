"""
Plot Semi-Major Axis vs. Eccentricity.
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import kepler_helpers as kh
import other_helpers as oh
import pandas_helpers as ph
import constants as C
import brewer2mpl as b2m
import argparse
import os
from time import gmtime, strftime
from glob import glob

# Load Colors
c3 = b2m.get_map('Dark2', 'Qualitative', 3)

# Scale Scatter
m, n = oh.mkline(5.0/2000.0, 1.0, 2.0, 36.0)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--tag", \
                    help="Title Tag. Defaults Current Working Dir.")
parser.add_argument("--cidx", type=int, default=0, \
                    help="Color Index (0, 1, or 2)")
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument("--all", action="store_true", \
                   help="Use All Snapshots.")
group1.add_argument("--custom", type=int, nargs="+", \
                   help="Use Custom Snapshot Range.")
args = parser.parse_args()

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Output set must be defined by three numbers."
        sys.exit()

# Full Set. Based On First/Count_In Directory.
if args.all:
    nsteps = []
    globs = glob("Out*.dat")
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)
    print "// Reading %i Outputs" % len(nsteps)

# Custom Set
if args.custom:
    # Build Output Number Array (From Input)
    nsteps = \
        np.mgrid[args.custom[0]:args.custom[1]+args.custom[2]:args.custom[2]]
    print "// Using Outputs %012d:%012d:%012d (%i Total)" % \
        ( args.custom[0], args.custom[1], args.custom[2], len(nsteps) )

# Set Unset Tag
if not args.tag:
    args.tag = os.getcwd()

# Determine Run Names
# Magic: Glob all ouputs, pick first, split off path,
#        extract filename, strip snap/#, strip "Out", et voila.
# Cf. http://www.webmasterwords.com/python-split-and-join-examples
run_names = []
globs = glob("Out*.dat")
globs = sorted(globs)
fname = globs[0].split("/")[-1]
run_name = "_".join(fname.split("_")[:-1])[3:]

# Loop Steps
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Load
    fname = "Out%s_%012d.dat" % (run_name, nstep)
    df = ph.read_output(fname)

    # Plot
    fig, ax = plt.subplots(1,1)
    s = df[df.mass<5.0].mass * m + n
    ax.scatter(df[df.mass<5.0].a, df[df.mass<5.0].e, \
               s=s**2.0, \
               color=c3.mpl_colors[args.cidx], \
               edgecolor=c3.mpl_colors[args.cidx])

    # Style
    ax.set_xlim([0,6])
    ax.set_ylim([0,1])
    ax.set_xlabel("Semi-Major Axis (AU)")
    ax.set_ylabel("Eccentricity")
    ax.set_title("%s" % args.tag)

    # Save Figure
    fig.savefig("ae_%012d.png" % nstep)

    # Clean Up
    plt.close(fig)
    del df

