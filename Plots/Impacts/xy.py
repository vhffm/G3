"""
Plot XY.
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
import sys
from time import gmtime, strftime
from glob import glob

# Load Colors
c3 = b2m.get_map('Dark2', 'Qualitative', 3)

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--no_orbits", type="store_true", \
                    help="Do Not Plot Ellipses.")
parser.add_argument("--count_in", type=int, \
                    help='Determine Output Range from This Line In Dirlist')
parser.add_argument("--tag", \
                    help="Title Tag. Defaults Current Working Dir.")
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument("--all", action="store_true", \
                   help="Use All Snapshots.")
group1.add_argument("--custom", type=int, nargs="+", \
                   help="Use Custom Snapshot Range.")
args = parser.parse_args()

# List of Directories.
# Format:
# run_tag_01,/path/to/director/01
# run_tag_02,/path/to/director/02
# ...
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    run_tags = []
    for line in lines:
        dirs.append(line.split(",")[1])
        run_tags.append(line.split(",")[0])
    print "// Reading %i Directories" % len(dirs)

# Sanity Check
if args.custom:
    if not len(args.custom) == 3:
        print "!! Output set must be defined by three numbers."
        sys.exit()

# Full Set. Based On First/Count_In Directory.
if args.all:
    nsteps = []
    if args.count_in:
        idir = args.count_in - 1
    else:
        idir = 0
    globs = glob("%s/Out*.dat" % dirs[idir])
    globs = sorted(globs)
    for g in globs:
        nstep = int(g.split(".")[0].split("_")[-1])
        nsteps.append(nstep)
    print "// Scanned %s For Output Range" % dirs[idir]
    print "// Reading %i Outputs Per Directory" % len(nsteps)

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
for idir, cdir in enumerate(dirs):
    globs = glob("%s/Out*.dat" % cdir)
    globs = sorted(globs)
    fname = globs[0].split("/")[-1]
    run_name = "_".join(fname.split("_")[:-1])[3:]
    run_names.append(run_name)

# Loop Steps
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Load
    fnames = []
    for idir, cdir in enumerate(dirs):
        fnames.append("%s/Out%s_%012d.dat" % \
                      (cdir, run_names[idir], nstep))
    df = ph.read_output_and_stack(fnames)

    # Ellipses
    if not args.no_orbits:
        xell, yell, _ = kh.compute_ellipseX(df[df.mass>0.0].a, \
                                            df[df.mass>0.0].e, \
                                            df[df.mass>0.0].i, \
                                            df[df.mass>0.0].Omega, \
                                            df[df.mass>0.0].omega)

    # Plot
    fig, ax = plt.subplots(1,1)

    # Hex Density
    ax.hexbin(df[df.mass==0.0].x, \
              df[df.mass==0.0].y, \
              cmap=mpl.cm.bone_r, \
              mincnt=1, \
              bins="log", \
              gridsize=512, \
              extent=[-60, 60, -60, 60], \
              edgecolors="grey", \
              linewidths=0.0, \
              alpha=1.0)

    # Plot Ellipses
    if not args.no_orbits:
        ax.plot(xell[:-1,:].T, \
                yell[:-1,:].T, \
                color="k", linewidth=1.0, alpha=0.5)

    # Style
    ax.set_aspect("equal")
    ax.set_xlim([-60,60])
    ax.set_ylim([-60,60])
    ax.set_xlabel("X (AU)")
    ax.set_ylabel("Y (AU)")

    # Titles
    # ax.set_title("%s" % args.tag)
    ax.set_title("%s *** %.2e yr *** %012d steps" % \
                 (args.tag, df.head(1).time, nstep))

    # Save Figure
    if args.no_orbits:
        fig.savefig("xyNoEllipses_%012d.png" % nstep)
    else:
        fig.savefig("xy_%012d.png" % nstep)

    # Clean Up
    plt.close(fig)
    del df
