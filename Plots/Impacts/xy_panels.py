"""
Plot XY. 6 Panels. Test and Massive Particles.

Plot Hexbin Particle Density for Test Particles
Plot Orbits for Massive Particles.

Speedup/slowdown wrt HDF5 reading varies. Probaly FS related.
Still save on storage space and pre-processing.
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import constants as C
import argparse
import sys
import os
import pandas as pd
from glob import glob
from time import gmtime, strftime
import kepler_helpers as kh

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument("--count_in", type=int, \
                    help='Determine Output Range from This Line In Dirlist')
parser.add_argument("--tag", \
                    help='Title Tag. Defaults Current Working Dir.')
group1 = parser.add_mutually_exclusive_group(required=True)
group1.add_argument('--all', action='store_true', \
                   help="Use All Snapshots.")
group1.add_argument('--custom', type=int, nargs='+', \
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

# Determine Run Name
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

# CSV Format
names_cols = [ "time", "pid", "mass", "radius", \
               "x", "y", "z", "vx", "vy", "vz", \
               "Sx", "Sy", "Sz", \
               "amin", "amax", "emin", "emax", \
               "aecount", "aecountT", "enccount", "test", "X" ]
touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9 ]
types_cols = { "time": np.float64, "pid": np.int64, "mass": np.float64, \
               "radius": np.float64, \
               "x": np.float64, "y": np.float64, "z": np.float64, \
               "vx": np.float64, "vy": np.float64, "vz": np.float64, \
               "Sx": np.float64, "Sy": np.float64, "Sz": np.float64, \
               "amin": np.float32, "amax": np.float32, \
               "emin": np.float32, "emax": np.float32, \
               "aecount": np.int64, "aecountT": np.int64, \
               "enccount": np.int64, "test": np.float64, "X": np.float32 }

#
# Loop Me, Baby
# https://www.youtube.com/watch?v=DVpq59F9I9o
# 
last_valid_step = [ 0, 0, 0, 0, 0, 0 ]
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(2,3)
    fig.set_size_inches(12,6)
    axflat = axarr.flatten()

    # Loop Runs
    for idir, cdir in enumerate(dirs):
        # Load Data
        try:
            fname = "%s/Out%s_%012d.dat" % (cdir, run_names[idir], nstep)
            df = pd.read_csv(fname, sep=" ", \
                             header=None, names=names_cols, dtype=types_cols, \
                             usecols=touse_cols, \
                             index_col=1)

            # Extract Time
            tout = df.head(1).time

            # Record Last Valid Step
            last_valid_step[idir] = nstep

            # Make Copies
            df0 = df[df.mass==0].copy()
            dfm = df[df.mass>0].copy()

            # Cast to Numpy
            x = np.asarray(dfm.x)
            y = np.asarray(dfm.y)
            z = np.asarray(dfm.z)
            vx = np.asarray(dfm.vx)
            vy = np.asarray(dfm.vy)
            vz = np.asarray(dfm.vz)
            m = np.asarray(dfm.mass)

            # Convert Frame
            x, vx = kh.helio2bary(x, vx, m)
            y, vy = kh.helio2bary(y, vy, m)
            z, vz = kh.helio2bary(z, vz, m)

            # Compute Keplerian Elements, Append Columns
            dfm["a"], dfm["e"], dfm["i"], dfm["Omega"], dfm["omega"], _ = \
                kh.cart2kepX(x, y, z, vx, vy, vz, m)

            # Convert Mass
            # dfm.mass *= C.msun/C.mearth
            
            # Plot Test Particles
            axflat[idir].hexbin(df0.x, \
                                df0.y, \
                                cmap=mpl.cm.bone_r, \
                                mincnt=1, \
                                bins="log", \
                                gridsize=128, \
                                extent=[-60, 60, -60, 60], \
                                edgecolors="grey", \
                                linewidths=0.0, \
                                alpha=1.0)
            
            # Compute Ellipses
            xell, yell, _ = \
                kh.compute_ellipseX(dfm.a, \
                                    dfm.e, \
                                    dfm.i, \
                                    dfm.Omega, \
                                    dfm.omega)
                
            # Plot Planets
            axflat[idir].plot(xell.T, yell.T, color="k", linewidth=1.0)
            axflat[idir].scatter(dfm.x, \
                                 dfm.y, \
                                 c="k", marker="+", s=10**2, \
                                 alpha=1.0, zorder=3)
            
            # Title
            axflat[idir].set_title("%s" % run_tags[idir])

            # Clean Up
            del df

        except:
            fname = "%s/Out%s_%012d.dat" % (cdir, run_names[idir], \
                                            last_valid_step[idir])
            df = pd.read_csv(fname, sep=" ", \
                             header=None, names=names_cols, dtype=types_cols, \
                             usecols=touse_cols, \
                             index_col=1)

            # Make Copies
            df0 = df[df.mass==0].copy()
            dfm = df[df.mass>0].copy()

            # Cast to Numpy
            x = np.asarray(dfm.x)
            y = np.asarray(dfm.y)
            z = np.asarray(dfm.z)
            vx = np.asarray(dfm.vx)
            vy = np.asarray(dfm.vy)
            vz = np.asarray(dfm.vz)
            m = np.asarray(dfm.mass)

            # Convert Frame
            x, vx = kh.helio2bary(x, vx, m)
            y, vy = kh.helio2bary(y, vy, m)
            z, vz = kh.helio2bary(z, vz, m)

            # Compute Keplerian Elements, Append Columns
            dfm["a"], dfm["e"], dfm["i"], dfm["Omega"], dfm["omega"], _ = \
                kh.cart2kepX(x, y, z, vx, vy, vz, m)

            # Convert Mass
            # dfm.mass *= C.msun/C.mearth
            
            # Plot Test Particles
            axflat[idir].hexbin(df0.x, \
                                df0.y, \
                                cmap=mpl.cm.bone_r, \
                                mincnt=1, \
                                bins="log", \
                                gridsize=128, \
                                extent=[-60, 60, -60, 60], \
                                edgecolors="grey", \
                                linewidths=0.0, \
                                alpha=1.0)
            
            # Compute Ellipses
            xell, yell, _ = \
                kh.compute_ellipseX(dfm.a, \
                                    dfm.e, \
                                    dfm.i, \
                                    dfm.Omega, \
                                    dfm.omega)
                
            # Plot Planets
            axflat[idir].plot(xell.T, yell.T, color="r", linewidth=1.0)
            axflat[idir].scatter(dfm.x, \
                                 dfm.y, \
                                 c="r", marker="+", s=10**2, \
                                 alpha=1.0, zorder=3)
            
            # Title
            axflat[idir].set_title("%s" % run_tags[idir])

            # Clean Up
            del df

        # Set Scale
        # for ax in axflat:
        #     ax.set_xscale("symlog", xlinthresh=0.1)
        #     ax.set_yscale("symlog", ylinthresh=0.1)

        # Set Limits
        for ax in axflat:
            ax.set_aspect("equal")
            ax.set_xlim([-60,60])
            ax.set_ylim([-60,60])
            
        # Remove Useless Labels
        for ii in [ 1, 2, 4, 5 ]:
            plt.setp(axflat[ii].get_yticklabels(), visible=False)

        for ii in [ 0, 1, 2 ]:
            plt.setp(axflat[ii].get_xticklabels(), visible=False)
            
        # Some Labels
        axarr[1,0].set_xlabel("X (AU)")
        axarr[1,0].set_ylabel("Y (AU)")
            
        # Fiture Title
        fig.suptitle("%s *** %.2e yr *** %012d steps" % \
            (os.getcwd(), tout, nstep))
    
    fig.savefig("xy_panels_%012d.png" % nstep)
    plt.close(fig)

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
