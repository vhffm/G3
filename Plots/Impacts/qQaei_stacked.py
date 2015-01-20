"""
Plot 3x3 Panels.
Horizontal: Three Bins in Semi-Major Axs at Time=0
Vertical  : qQ (Peri. vs Aphelion)
            ae (Semi-Major Axis vs. Eccentricity)
            qi (Perihelion vs. Inclination)
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import constants as C
import brewer2mpl as b2m
from glob import glob
import pandas as pd
import os
import sys
import logical_helpers as lh
import kepler_helpers as kh
import argparse
from time import gmtime, strftime
import logical_helpers as lh

# Global Styling
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['ytick.minor.size'] = 2
mpl.rcParams['xtick.labelsize']  = 'x-small'
mpl.rcParams['ytick.labelsize']  = 'x-small'
mpl.rcParams['axes.labelsize']   = 'x-small'
mpl.rcParams['axes.titlesize']   = 'small'
mpl.rcParams['font.size']        = 10

# Load Colors
c3 = b2m.get_map('Dark2', 'Qualitative', 3)

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
    print "// Using Outputs %012d:%012d:%012d" % \
        ( args.custom[0], args.custom[1], args.custom[2] )

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

# --------------
# Load Ejections
# --------------
print "// Loading Ejections"
names_cols = [ "time", \
               "index", "m", "r", "x", "y", "z", \
               "vx", "vy", "vz", "Sx", "Sy", "Sz", "case" ]
touse_cols = [ 0, 1, 13 ]
types_cols = { "index": np.int32, "case": np.int32 }

dfeject = pd.DataFrame()
for idir, cdir in enumerate(dirs):
    dfx = pd.read_csv("%s/Ejections%s.dat" % (cdir, run_names[idir]), \
                      sep=" ", \
                      header=None, names=names_cols, dtype=types_cols, \
                      usecols=touse_cols)
    #dfx.mass *= C.msun/C.mearth
    dfeject = dfeject.append(dfx, ignore_index=True)

# ---------------
# Load Collisions
# ---------------
print "// Loading Collisions"
names_cols = [ "time", \
               "indexi", "mi", "ri", "xi", "yi", "zi", "vxi", "vyi", "vzi", "Sxi", "Syi", "Szi",\
               "indexj", "mj", "rj", "xj", "yj", "zj", "vxj", "vyj", "vzj", "Sxj", "Syj", "Szj", "X" ]
# touse_cols = [ 0, 1, 2, 4, 5, 6, 13, 14, 16, 17, 18 ]
touse_cols = [ 0, 1, 2, 13, 14 ]
types_cols = { "indexi": np.int32, "indexj": np.int32 }

# Load data
dfcoll = pd.DataFrame()
for idir, cdir in enumerate(dirs):
    dfx = pd.read_csv("%s/Collisions%s.dat" % (cdir, run_names[idir]), \
                      sep=" ", \
                      header=None, names=names_cols, dtype=types_cols, \
                      usecols=touse_cols)
    dfx.mi *= C.msun/C.mearth
    dfx.mj *= C.msun/C.mearth
    dfcoll = dfcoll.append(dfx, ignore_index=True)

# ----------------
# Load Coordinates
# ----------------
names_cols = [ "time", "pid", "mass", "radius", "x", "y", "z", "vx", "vy", "vz", "Sx", "Sy", "Sz", \
               "amin", "amax", "emin", "emax", "aecount", "aecountT", "enccount", "test", "X" ]
touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9 ]
types_cols = { "pid": np.int32 }

#
# Loop me
#
diverged = False
for istep, nstep in enumerate(nsteps):
    print "// (%s UTC) Processing Snapshot %012d/%012d" % \
        (strftime("%H:%M:%S", gmtime()), nstep, nsteps[-1])

    # Loop Runs
    df = pd.DataFrame()
    for idir, cdir in enumerate(dirs):
        dfx = pd.read_csv("%s/Out%s_%012d.dat" % \
                          (cdir, run_names[idir], nstep), \
                          sep=" ", \
                          header=None, names=names_cols, dtype=types_cols, \
                          usecols=touse_cols, \
                          index_col=1)
        #dfx.mass *= C.msun/C.mearth
        #df = df.append(dfx, ignore_index=True)
        df = df.append(dfx)

    # Drop duplicates
    df.drop_duplicates(inplace=True)

    # Check if massive planets are (still) in sync
    if not diverged:
        if np.sum(df.drop_duplicates().index.isin([0,1,2,3,4,5,6,7,8])) > 9:
            diverged = True

    # Drop duplicates by index
    # http://stackoverflow.com/questions/13035764/remove-rows-with-duplicate-indices-pandas-dataframe-and-timeseries
    df["index"] = df.index
    df.drop_duplicates(subset=["index"], inplace=True)
    del df["index"]
    # df.sort(inplace=True)

    # Use np.asarray() to extract numpy arrays from the dataframe
    x = np.asarray(df.x); y = np.asarray(df.y); z = np.asarray(df.z)
    vx = np.asarray(df.vx); vy = np.asarray(df.vy); vz = np.asarray(df.vz)
    m = np.asarray(df.mass)

    # Convert to barycentric coordinates
    x, vx = kh.helio2bary(x, vx, m)
    y, vy = kh.helio2bary(y, vy, m)
    z, vz = kh.helio2bary(z, vz, m)

    # Compute Keplerian elements, and append the columns to the dataframe 
    df["a"], df["e"], df["i"], df["Omega"], df["omega"], df["M"] = \
        kh.cart2kepX(x, y, z, vx, vy, vz, m)

    # Compute Peri- and Aphelia
    df["Q"] = df["a"] * ( 1 + df["e"] )
    df["q"] = df["a"] * ( 1 - df["e"] )

    # Convert mass
    # df.mass *= C.msun/C.mearth

    # Output Keplerian elements for the first 15 particles
    # print len(df)
    # print df[["mass", "a", "e", "i", "Q", "q"]].head(15)

    # Compute Semi-Major Axis Bin Limits
    # NB: Assuming @IC Particle ID Increases with Semi-Major Axis
    if istep == 0:
        if np.sum(df.index.isin([25000, 80000])) < 2:
            print "!! Bin Boundary Particles Missing. Setting NaN."
            abins_lo = np.ones(3) * np.nan
            abins_hi = np.ones(3) * np.nan
        else:
            abins_lo = np.asarray([ \
                df[df.index>9].sort(columns=["a"]).head(1).iloc[0].loc["a"], \
                df[df.index==25000].iloc[0].loc["a"], \
                df[df.index==80000].iloc[0].loc["a"] \
            ])
            abins_hi = np.asarray([ \
                df[df.index==25000].iloc[0].loc["a"], \
                df[df.index==80000].iloc[0].loc["a"], \
                df[df.index>9].sort(columns=["a"]).tail(1).iloc[0].loc["a"] \
            ])

    # Initialize Filters
    # Centr Covers Jupiter/Saturn
    bool_inner = lh.xlogical_and([df.mass==0, df.index<25000, \
                                  df.e>=0.0, df.e<1.0])
    bool_centr = lh.xlogical_and([df.mass==0, df.index>=25000, \
                                  df.index<80000, df.e>=0.0, df.e<1.0])
    bool_outer = lh.xlogical_and([df.mass==0, df.index>=80000, \
                                  df.e>=0.0, df.e<1.0])
    # Startup-speak ("filter" is a reserved keyword)
    filtr = [ bool_inner, bool_centr, bool_outer ]

    # Plot
    fig, axmat = plt.subplots(3,3)
    fig.set_size_inches(12,9)

    # ------------------
    # Peri- vs. Aphelion
    # ------------------
    axarr = axmat[0]

    # Loop Axis
    for iax, ax in enumerate(axarr):

        # Inner Solar System
        ax.fill_between([1.0e-1, 1.0e6], [1.0e-1, 1.0e-1], [1.7, 1.7], \
                        facecolor=c3.mpl_colors[1], alpha=0.05, lw=0.5)
        ax.fill_between([1.0e-1, 1.0e6], [1.0e-1, 1.0e-1], [1.1, 1.1], \
                        facecolor=c3.mpl_colors[1], alpha=0.05, lw=0.5)
        ax.fill_between([1.0e-1, 1.7], [1.0e-1, 1.0e-1], [1.0e2, 1.0e2], \
                        facecolor=c3.mpl_colors[1], alpha=0.05, lw=0.5)
        ax.fill_between([1.0e-1, 1.1], [1.0e-1, 1.0e-1], [1.0e2, 1.0e2], \
                        facecolor=c3.mpl_colors[1], alpha=0.05, lw=0.5)

        # Massive
        ax.scatter(df.Q[df.mass>0], df[df.mass>0].q, \
                   s=(df.mass[df.mass>0]/(C.mmercury/C.msun)+30)**(2./3.), \
                   c=c3.mpl_colors[1], alpha=0.8, lw=0.5)

        # Test Particles
        ax.hexbin(df[filtr[iax]].Q, df[filtr[iax]].q, \
                  cmap=mpl.cm.bone_r, \
                  mincnt=1, \
                  bins="log", \
                  gridsize=256, \
                  xscale="log", \
                  yscale="log", \
                  extent=[-1, 6, -1, 2], \
                  edgecolors="grey", \
                  linewidths=0.0, \
                  alpha=1.0, \
                  vmin=-0.5)

        # Reference Circular Orbit Line
        ax.plot([0.1,1,10,100], [0.1,1,10,100], c='k', alpha=0.2, lw=0.5)

        # Style
        ax.set_xlim([1.0e-1,1.0e6])
        ax.set_ylim([1.0e-1,1.0e2])
        ax.set_xlabel("Q=a(1+e) (AU)")
        if iax == 0:
            ax.set_ylabel("q=a(1-e) (AU)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.xaxis.labelpad = 0.5
        
        # Test Particles Counts
        ax.text(0.5, 0.38, r"$N = %i$" % \
                np.sum(filtr[iax]), \
                horizontalalignment='left', color='black', \
                transform=ax.transAxes)
        ax.text(0.5, 0.30, r"$N(q<1.1,i<5) = %i$" % \
                np.sum(np.logical_and(df[filtr[iax]].q<1.1, df[filtr[iax]].i<5.0*C.d2r)), \
                horizontalalignment='left', color='black', \
                transform=ax.transAxes)
        ax.text(0.5, 0.22, r"$N(q<1.7,i<5) = %i$" % \
                np.sum(np.logical_and(df[filtr[iax]].q<1.7, df[filtr[iax]].i<5.0*C.d2r)), \
                horizontalalignment='left', color='black', \
                transform=ax.transAxes)
        ax.text(0.5, 0.14, r"$N(q<1.7) = %i$" % \
                np.sum(df[filtr[iax]].q<1.7), \
                horizontalalignment='left', color='black', \
                transform=ax.transAxes)
        ax.text(0.5, 0.06, r"$N(q<1.1) = %i$" % \
                np.sum(df[filtr[iax]].q<1.1), \
                horizontalalignment='left', color='black', \
                transform=ax.transAxes)

        # Inner Planet Collisions
        if iax == 0:
            ax.text(0.02, 0.92, r"$N_\mathrm{Impacts}(M,V,E,M) = (%i,%i,%i,%i)$" % \
                    (np.sum(np.logical_and(dfcoll.time<df.iloc[0].loc["time"], np.logical_or(dfcoll.indexi==0, dfcoll.indexj==0))), \
                     np.sum(np.logical_and(dfcoll.time<df.iloc[0].loc["time"], np.logical_or(dfcoll.indexi==1, dfcoll.indexj==1))), \
                     np.sum(np.logical_and(dfcoll.time<df.iloc[0].loc["time"], np.logical_or(dfcoll.indexi==2, dfcoll.indexj==2))), \
                     np.sum(np.logical_and(dfcoll.time<df.iloc[0].loc["time"], np.logical_or(dfcoll.indexi==3, dfcoll.indexj==3)))), \
                    horizontalalignment='left', color='black', \
                    transform=ax.transAxes)

        # Divergence Warning
        if diverged:
            if iax > 0:
                ax.text(0.5, 0.92, "!!! DIVERGED !!!", \
                horizontalalignment='left', color='red', \
                transform=ax.transAxes)
                
        # Titles
        ax.set_title("%.2f < a(t=0) < %.2f AU" % (abins_lo[iax], abins_hi[iax]))

    # --------------------------------
    # Semi-Major Axis vs. Eccentricity
    # --------------------------------
    axarr = axmat[1]

    # Loop Axes
    for iax, ax in enumerate(axarr):
        
        # Inner Solar System
        ax.fill_between([1.0e-1, 1.0], [0, 0], [1, 1], \
                        facecolor=c3.mpl_colors[1], alpha=0.075, lw=0.5)
        ax.fill_between([1.0e-1, 1.5], [0, 0], [1, 1], \
                        facecolor=c3.mpl_colors[1], alpha=0.075, lw=0.5)
        
        # Massive
        ax.scatter(df.a[df.mass>0], df[df.mass>0].e, \
                   s=(df.mass[df.mass>0]/(C.mmercury/C.msun)+30)**(2./3.), \
                   c=c3.mpl_colors[1], alpha=0.8, lw=0.5)

        # Test Particles
        ax.hexbin(df[filtr[iax]].a, df[filtr[iax]].e, \
                  cmap=mpl.cm.bone_r, \
                  mincnt=1, \
                  bins="log", \
                  gridsize=256, \
                  xscale="log", \
        #           yscale="log", \
                  extent=[-1, 3, 0, 1], \
                  edgecolors="grey", \
                  linewidths=0.0, \
                  alpha=1.0, \
                  vmin=-0.5)

        # Style
        ax.set_xlim([1.0e-1,1.0e3])
        ax.set_ylim([0,1])
        ax.set_xlabel("a (AU)")
        if iax == 0:
            ax.set_ylabel("e")
        ax.xaxis.labelpad = 0.5
            
        # Collision/Ejection Counts
        if iax == 0:
            ax.text(0.02, 0.92, r"$N_{\mathrm{Total}} = %i$" % \
                    np.sum(df.mass==0), \
                    horizontalalignment='left', color='black', \
                    transform=ax.transAxes)
            ax.text(0.02, 0.84, r"$N_{\mathrm{Ejected}} = %i$" % \
                    np.sum(np.logical_and(dfeject.time<df.iloc[0].loc["time"], dfeject.case==-3)), \
                    horizontalalignment='left', color='black', \
                    transform=ax.transAxes)
            ax.text(0.02, 0.76, r"$N_{\mathrm{Infall}} = %i$" % \
                    np.sum(np.logical_and(dfeject.time<df.iloc[0].loc["time"], dfeject.case==-2)), \
                    horizontalalignment='left', color='black', \
                    transform=ax.transAxes)
            ax.text(0.02, 0.68, r"$N_{\mathrm{Collision}} = %i$" % \
                    np.sum(dfeject.time<df.iloc[0].loc["time"]), \
                    horizontalalignment='left', color='black', \
                    transform=ax.transAxes)

        # Divergence Warning
        if diverged:
            ax.text(0.5, 0.92, "!!! DIVERGED !!!", \
            horizontalalignment='left', color='red', \
            transform=ax.transAxes)

    # --------------------------
    # Perihelion vs. Inclination
    # --------------------------
    axarr = axmat[2]

    # Loop Axis
    for iax, ax in enumerate(axarr):
        
        # Inner Solar System
        ax.fill_between([1.0e-1, 1.0], [1.0e-2, 1.0e-2], [1.0e2, 1.0e2], \
                        facecolor=c3.mpl_colors[1], alpha=0.075, lw=0.5)
        ax.fill_between([1.0e-1, 1.5], [1.0e-2, 1.0e-2], [1.0e2, 1.0e2], \
                        facecolor=c3.mpl_colors[1], alpha=0.075, lw=0.5)
        
        # Massive
        ax.scatter(df.a[df.mass>0], df[df.mass>0].i * C.r2d, \
                   s=(df.mass[df.mass>0]/(C.mmercury/C.msun)+30)**(2./3.), \
                   c=c3.mpl_colors[1], alpha=0.8, lw=0.5)

        # Test Particles
        ax.hexbin(df[filtr[iax]].a, df[filtr[iax]].i * C.r2d, \
                  cmap=mpl.cm.bone_r, \
                  mincnt=1, \
                  bins="log", \
                  gridsize=256, \
                  xscale="log", \
                  yscale="log", \
                  extent=[-1, 3, -2, 2], \
                  edgecolors="grey", \
                  linewidths=0.0, \
                  alpha=1.0, \
                  vmin=-0.5)

        # Style
        ax.set_xlim([1.0e-1,1.0e3])
        ax.set_ylim([1.0e-2,1.0e2])
        ax.set_xlabel("a (AU)")
        if iax == 0:
            ax.set_ylabel("i (Degree)")
        ax.xaxis.labelpad = 0.5

        # Divergence Warning
        if diverged:
            ax.text(0.5, 0.92, "!!! DIVERGED !!!", \
            horizontalalignment='left', color='red', \
            transform=ax.transAxes)

    # Set Title
    fig.suptitle("%s *** %.2e yr *** %012d steps" % \
                 (os.getcwd(), df.head(1).time, nstep))
    fig.subplots_adjust(top=0.93)

    # Save Figure
    fig.savefig("qQaei_%012d.png" % nstep)
    plt.close(fig)

print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
