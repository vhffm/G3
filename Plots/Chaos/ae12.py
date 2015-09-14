"""
Plot Semi-Major Axis vs. Eccentricity. 12 Panels.
Requires Coordinates_XX.hdf5 and Resonances_XX.hdf5.

Dirlist Format:
/path/Coordinates_01.hdf5,/path/Resonances_01.hdf5
/path/Coordinates_02.hdf5,/path/Resonances_02.hdf5
...
/path/Coordinates_NN.hdf5,/path/Resonances_NN.hdf5
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import other_helpers as oh


# List of Directories
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    crd_files = []
    res_files = []
    for line in lines:
        line = line.strip().split(',')
        crd_files.append(line[0])
        res_files.append(line[1])
    print "// Reading %i Runs" % len(crd_files)

# Load Resonances
dfr_all = []
for res_file in res_files:
    with pd.HDFStore("%s" % res_file, 'r') as store:
        dfr = store['df']
        dfr_all.append(dfr)
        del dfr

# Load Coordinates
dfo_all = []
for crd_file in crd_files:
    with pd.HDFStore("%s" % (crd_file), 'r') as store:
        dfo = store['df']
        dfo_all.append(dfo)
        del dfo

# Determine Step Range
nsteps = np.asarray(dfo_all[0].nstep.unique())

# Loop Steps
print "// Processing %i Outputs per Run" % len(nsteps)
for nstep in nsteps:
    print "** Step %012d/%012d" % (nstep, nsteps[-1])

    # Setup Figure
    fig, axarr = plt.subplots(4,3)
    fig.set_size_inches(16,12)

    # Loop Runs
    for irow in [ 0, 1, 2, 3 ]:
        for icol in [ 0, 1, 2 ]:

            # Index Juggling
            irun = irow * 3 + icol

            # Pick Axis
            ax = axarr[irow,icol]

            if irun < len(dfo_all):

                # Pick Run
                dfo = dfo_all[irun]
                dfr = dfr_all[irun]

                # Pick Step
                dfo = dfo[dfo.nstep == nstep]
                dfr = dfr[dfr.nstep == nstep]

                # Extract Jupiter Semi-Major Axis
                a_j = dfo[dfo.pid == 2000].a.iloc[0]

                # Extract Jupiter Crossing Orbits
                Q_cross = a_j
                e_cross = np.linspace(0,1,128)
                a_cross = Q_cross / ( 1.0 + e_cross )

                # Remove Giant Planets
                dfo = dfo[dfo.mass < 5.0]

                # Scale Line
                ms, ns = oh.mkline(5.0/2000.0, 2.0, 2.0, 36.0)
                ss = np.asarray(dfo.mass * ms + ns)

                # #############################################
                # #############################################
                # PLANETESIMALS
                # #############################################
                # #############################################
                
                # Planetesimals
                ax.scatter(dfo.a, dfo.e, s=ss**2.0, \
                           c='k', edgecolor='none', alpha=0.8, zorder=0)

                # #############################################
                # #############################################
                # RESONANCE LOCATIONS
                # #############################################
                # #############################################

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
    for ax in axarr[:3,:].flatten():
        plt.setp(ax.get_xticklabels(), visible=False)
    for ax in axarr[:,1:].flatten():
        plt.setp(ax.get_yticklabels(), visible=False)

    # Axis Labels
    # for ax in axarr[-1,:].flatten():
    #     ax.set_xlabel('Semi-Major Axis (AU)')
    # for ax in axarr[:,0].flatten():
    #     ax.set_ylabel('Eccentricity')
    axarr[-1,0].set_xlabel('Semi-Major Axis (AU)')
    axarr[-1,0].set_ylabel('Eccentricity')
            
    # Title, Legend
    fig.suptitle("nstep = %.0e, time = %.2e Myr" % \
        (nstep, dfo.time.iloc[0]/1.0e6))

    # Join Panels
    fig.subplots_adjust(hspace=0,wspace=0,top=0.96)

    # Save Figure
    fig.savefig("ae12_%012d.png" % nstep)

    # Clean Up
    del ax
    del fig

