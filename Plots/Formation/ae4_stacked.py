"""
Plot Semi-Major Axis vs. Eccentricity. 4 Panels. Stack.
Requires Coordinates_XX.hdf5.

Dirlist Format:
tag_01,/path/Coordinates_01.hdf5
tag_01,/path/Coordinates_02.hdf5
tag_02,/path/Coordinates_01.hdf5
tag_02,/path/Coordinates_02.hdf5
...
tag_NN,/path/Coordinates_MM.hdf5
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import other_helpers as oh
import constants as C


# List of Directories
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    tags = []
    crd_files = []
    for line in lines:
        line = line.strip().split(',')
        tags.append(line[0])
        crd_files.append(line[1])
    print "// Reading %i Runs" % len(crd_files)

# Create Numerical IDs for Unique Tags
tag_ids = []; ktag = 0; first = True
for tag in tags:
    if not first:
        if not tag == tag_last:
            ktag += 1
    else:
        first = False
    tag_ids.append(ktag)
    tag_last = tag

# Retain Unique Tags
tags_unique = []
for tag in tags:
    if not tag in tags_unique:
        tags_unique.append(tag)

# Create Lists with Filenames to Load
fnames_all = []
fnames_loc = []
for itag_id, tag_id in enumerate(tag_ids):
    # Append filename to local
    fnames_loc.append(crd_files[itag_id])
    # Is the next tag the same as me or are we last?
    # If not, append to global and reset local. 
    # Otherwise, do nothting.
    if (itag_id == len(tag_ids)-1) or (not (tag_ids[itag_id+1] == tag_id)):
        fnames_all.append(fnames_loc)
        fnames_loc = []

# Load Coordinates
df_all = []
print "// Loading Coordinates Files"
for fnames_loc in fnames_all:
    df = pd.DataFrame()
    for fname in fnames_loc:
        print "** %s" % fname
        with pd.HDFStore("%s" % fname, 'r') as store:
            df_tmp = store['df']
            df = pd.concat([df, df_tmp])
            del df_tmp
    df_all.append(df.reset_index(drop=True))

# Determine Step Range
nsteps = np.asarray(df_all[0].nstep.unique())

# Global Mass Scaling
mx, nx = oh.mkline(5.0/8192.0, 1.0, 1.6, 12.0)
mc, nc = oh.mkline(5.0/8192.0, 0.6, 1.6, 1.0)

# Loop
print "// Processing %i Steps" % len(nsteps)
for nstep in nsteps:
    print "** Step %012d/%012d" % (nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(2,2)
    fig.set_size_inches(8,8)

    # Loop Sims
    for isim, df in enumerate(df_all):
        
        # Setup Axis
        ax = axarr.flatten()[isim]

        # Filter Step, Extract Time
        df = df[df.nstep == nstep]
        if isim == 0:
            tout = df.time.iloc[0]
        
        # Scaling
        s = mx * np.asarray(df.mass) + nx
        c = mc * np.asarray(df.mass) + nc

        # Plot
        ax.scatter(df.a, df.e, \
                   c=c, s=s**2, alpha=0.8, \
                   edgecolor='none', \
                   cmap=mpl.cm.Greys, vmin=0, vmax=1)
        
        # Style
        ax.set_xlim([0,6])
        ax.set_ylim([0,0.6])

        # Annotation
        ax.text(0.1, 0.9, "(%i, %i, %i)" % ( np.sum(df.mass>=C.m_cutoff), \
                                             np.sum(df.mass<C.m_cutoff), \
                                             len(df) ), \
                va='top', ha='left', transform=ax.transAxes, \
                fontsize='small')
        
        # Title
        ax.set_title("%s" % tags_unique[isim])
        
    # Remove Labels
    for iax, ax in enumerate(axarr.flatten()):
        if iax < 2:
            plt.setp(ax.get_xticklabels(), visible=False)
        if iax == 1 or iax == 3:
            plt.setp(ax.get_yticklabels(), visible=False)    
            
    # Axis Labels
    axarr.flatten()[2].set_xlabel('Semi-Major Axis (AU)')
    axarr.flatten()[2].set_ylabel('Eccentricity')

    # Main Title
    txtsup = "nstep = %012d // time = %.2e Myr" % (nstep, tout/1.0e6)
    txtsup += "  // (Embryos/Planets, Planetesimals, Total Particles)"
    fig.suptitle(txtsup)

    # Adjust Spacing
    # fig.subplots_adjust(top=0.82)

    # Save
    fig.savefig("ae4_stacked_%012d.png" % nstep, bbox_inches='tight')

    # Clean Up
    plt.close(fig)
