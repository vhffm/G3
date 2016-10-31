"""
Plot Semi-Major Axis vs. Eccentricity. 2 Panels. Stack. Wide.
Requires Coordinates_XX.hdf5, Resonances_XX.hdf5.

Dirlist Format:
tag_01,/path/Coordinates_01.hdf5,/path/Resonances_01.hdf5
tag_01,/path/Coordinates_02.hdf5,/path/Resonances_02.hdf5
tag_02,/path/Coordinates_01.hdf5,/path/Resonances_01.hdf5
tag_02,/path/Coordinates_02.hdf5,/path/Resonances_02.hdf5
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
    res_files = []
    for line in lines:
        line = line.strip().split(',')
        tags.append(line[0])
        crd_files.append(line[1])
        res_files.append(line[2])
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

# Create Lists with Filenames to Load -- Coordinates
fnames_crd_all = []
fnames_crd_loc = []
for itag_id, tag_id in enumerate(tag_ids):
    # Append filename to local
    fnames_crd_loc.append(crd_files[itag_id])
    # Is the next tag the same as me or are we last?
    # If not, append to global and reset local. 
    # Otherwise, do nothting.
    if (itag_id == len(tag_ids)-1) or (not (tag_ids[itag_id+1] == tag_id)):
        fnames_crd_all.append(fnames_crd_loc)
        fnames_crd_loc = []

# Create Lists with Filenames to Load -- Resonances
fnames_res_all = []
fnames_res_loc = []
for itag_id, tag_id in enumerate(tag_ids):
    # Append filename to local
    fnames_res_loc.append(res_files[itag_id])
    # Is the next tag the same as me or are we last?
    # If not, append to global and reset local. 
    # Otherwise, do nothting.
    if (itag_id == len(tag_ids)-1) or (not (tag_ids[itag_id+1] == tag_id)):
        fnames_res_all.append(fnames_res_loc)
        fnames_res_loc = []

# Load Coordinates
dfo_all = []
print "// Loading Coordinates Files"
for fnames_crd_loc in fnames_crd_all:
    dfo = pd.DataFrame()
    for fname in fnames_crd_loc:
        print "** %s" % fname
        with pd.HDFStore("%s" % fname, 'r') as store:
            dfo_tmp = store['df']
            dfo = pd.concat([dfo, dfo_tmp])
            del dfo_tmp
    dfo_all.append(dfo.reset_index(drop=True))

# Load Resonances
dfr_all = []
print "// Loading Resonance Files"
for fnames_res_loc in fnames_res_all:
    dfr = pd.DataFrame()
    if not fnames_res_loc[0] == 'NONE':
        for fname in fnames_res_loc:
            print "** %s" % fname
            with pd.HDFStore("%s" % fname, 'r') as store:
                dfr_tmp = store['df']
                dfr = pd.concat([dfr, dfr_tmp])
                del dfr_tmp
        # Compute Mean Location for Resonances, MultiIndex Magic
        dfr.set_index(['nrun', 'nstep'], drop=True, append=False, inplace=True)
        dfr = dfr.mean(level='nstep')
        dfr['nstep'] = dfr.index
        dfr.reset_index(drop=True, inplace=True)
    dfr_all.append(dfr)

# Determine Step Range
nsteps = np.asarray(dfo_all[0].nstep.unique())

# Global Mass Scaling
ms, ns = oh.mkline(5.0/2000.0, 1.0, 2.0, 16.0)

# Loop Steps
print "// Processing %i Outputs per Run" % len(nsteps)
for nstep in nsteps:
    print "** Step %012d/%012d" % (nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(1,2)
    fig.set_size_inches(12,3)

    # Loop Runs
    for isim in [ 0, 1 ]:

        # Select Axis
        ax = axarr[isim]

        # Select Sim
        dfo = dfo_all[isim]
        dfr = dfr_all[isim]

        # Pick Step
        dfo = dfo[dfo.nstep == nstep]
        if not fnames_res_all[isim][0] == 'NONE':
            dfr = dfr[dfr.nstep == nstep]

        # Extract Jupiter (Mean) Semi-Major Axis
        if np.sum(np.asarray(dfo.pid == 2000)) > 0:
            a_j = np.mean(np.asarray(dfo[dfo.pid == 2000].a))
        else:
            a_j = np.nan
        
        # Extract Jupiter Crossing Orbits
        Q_cross = a_j
        e_cross = np.linspace(0,1,128)
        a_cross = Q_cross / ( 1.0 + e_cross )

        # Remove Giant Planets
        dfo = dfo[dfo.mass < 5.0]

        # Scaling
        s = ms * np.asarray(dfo.mass) + ns
        c = np.ones([len(dfo.mass),4]) * np.array([0.6,0.6,0.6,1.0])

        # #############################################
        # #############################################
        # PLANETESIMALS
        # #############################################
        # #############################################
        
        # Plot
        ax.scatter(dfo.a, dfo.e, \
                   c=c, s=s**2, alpha=0.7, \
                   edgecolor='k', \
                   linewidth=0.5)

        # #############################################
        # #############################################
        # RESONANCE LOCATIONS
        # #############################################
        # #############################################

        if not fnames_res_all[isim][0] == 'NONE':

            # Resonances (Secular)
            ax.plot([dfr.nu_5.iloc[0], dfr.nu_5.iloc[0]], \
                    [0, 1], \
                    'k', lw=0.5, label='nu_5', zorder=-10)
            ax.plot([dfr.nu_6.iloc[0], dfr.nu_6.iloc[0]], \
                    [0, 1], \
                    'k', lw=0.5, label='nu_6', zorder=-10)
            if not np.isnan(dfr.nu_16_01.iloc[0]):
                ax.plot([dfr.nu_16_01.iloc[0], dfr.nu_16_01.iloc[0]], \
                        [0, 1], \
                        'k', lw=0.5, label='nu_16', zorder=-10)
            if not np.isnan(dfr.nu_16_02.iloc[0]):
                ax.plot([dfr.nu_16_02.iloc[0], dfr.nu_16_02.iloc[0]], \
                        [0, 1], \
                        'k', lw=0.5, label='nu_16', zorder=-10)

            # Resonances (MMR)
            ax.plot([dfr.two_to_one.iloc[0], dfr.two_to_one.iloc[0]], \
                    [0, 1], \
                    color=(0.6,0.6,0.6), lw=0.5, label='2:1', zorder=-10)
            ax.plot([dfr.three_to_one.iloc[0], dfr.three_to_one.iloc[0]], \
                    [0, 1], \
                    color=(0.6,0.6,0.6), lw=0.5, label='3:1', zorder=-10)
            ax.plot([dfr.three_to_two.iloc[0], dfr.three_to_two.iloc[0]], \
                    [0, 1], \
                    color=(0.6,0.6,0.6), lw=0.5, label='3:2', zorder=-10)

        # #############################################
        # #############################################
        # ANNOTATIONS
        # #############################################
        # #############################################

        if not fnames_res_all[isim][0] == 'NONE':

            # Annotate (Secular)
            ax.text(dfr.nu_5.iloc[0], 0.55, \
                    r'$\nu_5$', transform=ax.transData, \
                    va='center', ha='center', \
                    bbox=dict(facecolor='w', lw=0.5), fontsize='small')

            ax.text(dfr.nu_6.iloc[0], 0.55, \
                    r'$\nu_6$', transform=ax.transData, \
                    va='center', ha='center', \
                    bbox=dict(facecolor='w', lw=0.5), fontsize='small')
            
            if not np.isnan(dfr.nu_16_01.iloc[0]):
                ax.text(dfr.nu_16_01.iloc[0], 0.55, \
                        r'$\nu_{16}$', transform=ax.transData, \
                        va='center', ha='center', \
                        bbox=dict(facecolor='w', lw=0.5), fontsize='small')
            if not np.isnan(dfr.nu_16_02.iloc[0]):
                ax.text(dfr.nu_16_02.iloc[0], 0.55, \
                        r'$\nu_{16}$', transform=ax.transData, \
                        va='center', ha='center', \
                        bbox=dict(facecolor='w', lw=0.5), fontsize='small')

            # Annotate (MMR)
            ax.text(dfr.two_to_one.iloc[0], 0.475, \
                    r'$2:1$', transform=ax.transData, \
                    va='center', ha='center', \
                    bbox=dict(facecolor='w', lw=0.5), fontsize='x-small')

            ax.text(dfr.three_to_one.iloc[0], 0.475, \
                    r'$3:1$', transform=ax.transData, \
                    va='center', ha='center', \
                    bbox=dict(facecolor='w', lw=0.5), fontsize='x-small')

            ax.text(dfr.three_to_two.iloc[0], 0.475, \
                    r'$3:2$', transform=ax.transData, \
                    va='center', ha='center', \
                    bbox=dict(facecolor='w', lw=0.5), fontsize='x-small')

        # Number of Particles
        ax.text(0.1, 0.7, "(%i, %i, %i)" % ( np.sum(dfo.mass>=C.m_cutoff), \
                                             np.sum(dfo.mass<C.m_cutoff), \
                                             len(dfo) ), \
                va='top', ha='left', transform=ax.transAxes, \
                fontsize='small')

        # #############################################
        # #############################################
        # CUTOFF, JUPITER CROSSING
        # #############################################
        # #############################################

        # Jupiter Crossing
        # ax.plot(a_cross, e_cross, 'k--')
        ax.fill_between(a_cross, e_cross, 0.6, \
                        lw=0, facecolor='r', \
                        alpha=0.1, zorder=-20)

        # Inner Cutoff
        ax.fill_between([0.0, 0.1], 0.0, 0.6, \
                        lw=0, facecolor='r', \
                        alpha=0.1, zorder=20)

        # Title
        ax.set_title("%s" % tags_unique[isim], fontsize='small')

    # #############################################
    # #############################################
    # STYLING
    # #############################################
    # #############################################
            
    # Limits
    for ax in axarr.flatten():
        ax.set_xlim([0,5])
        ax.set_ylim([0,0.6])
        
    # Ticks
    for ax in axarr.flatten():
        ax.xaxis.set_major_locator(\
            mpl.ticker.MaxNLocator(prune='both', nbins='5'))
        ax.yaxis.set_major_locator(\
            mpl.ticker.MaxNLocator(prune='both', nbins='5'))
    
    # Hide Tick Labels
    # for ax in axarr[:3,:].flatten():
    #     plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axarr[1:].flatten():
        plt.setp(ax.get_yticklabels(), visible=False)

    # Axis Labels
    # for ax in axarr[-1,:].flatten():
    #     ax.set_xlabel('Semi-Major Axis (AU)')
    # for ax in axarr[:,0].flatten():
    #     ax.set_ylabel('Eccentricity')
    for ax in axarr:
        ax.set_xlabel('Semi-Major Axis (AU)', fontsize='small')
    axarr[0].set_ylabel('Eccentricity', fontsize='small')
            
    # Title, Legend
    tout = dfo.time.iloc[0]
    txtsup = "nstep = %012d // time = %.2e Yr" % (nstep, tout)
    txtsup += "  // (Embryos/Planets, Planetesimals, Total Particles)"
    fig.suptitle(txtsup, fontsize='small')

    # Join Panels
    # fig.subplots_adjust(hspace=0,wspace=0,top=0.96)
    fig.subplots_adjust(top=0.84)

    # Save Figure
    fig.savefig("ae2_stacked_wide_%012d.png" % nstep, bbox_inches='tight')

    # Clean Up
    plt.close(fig)
