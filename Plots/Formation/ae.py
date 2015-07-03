"""
Plot Semi-Major Axis vs. Eccentricity. Scale Point Size by Mass.
"""

import matplotlib as mpl; mpl.use('agg')
import matplotlib.pyplot as plt
import io_helpers as ioh
import glob
import argparse
import other_helpers as oh


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--run_name', default='run_09', \
                    help='Name of Simulation Run.')
parser.add_argument('--dir_name', default='.', \
                    help='Directory')
args = parser.parse_args()

# Glob Outputs
nsteps = []
globs = glob.glob("%s/Out*.dat" % args.dir_name)
globs = sorted(globs)
for g in globs:
    nstep = int(g.strip()[:-4].split("_")[-1])
    nsteps.append(nstep)

# #############################################################################
# Find Smallest/Largest Mass
# #############################################################################

# 0 - Smallest Mass from First Output
df = ioh.read_output("%s/Out_%s_%012d.dat" % (args.dir_name, \
                                              args.run_name, \
                                              nsteps[-1]), \
                     frame="heliocentric")
mmin = df["mass"].min()

# # 1 - Largest Mass from Last Output
# mmax = np.zeros(2)
# df = ioh.read_output("Out_%s_%012d.dat" % (args.run_name, nsteps[-1]), \
#                      frame="heliocentric")
# df = df[df.mass<16.0] 
# mmax[0] = df["mass"].max()
# del df

# # 1 - Largest Ejected Mass
# df = ioh.read_ejections_and_stack(["Ejections_%s.dat" % (args.run_name)])
# mmax[1] = df["m"].max()
# del df

# Set Scaling (Smallest/Largest)
# m, n = oh.mkline(mmin, 1.0, np.nanmax(mmax), 36.0)

# Set Scaling (Smallest/One-Earth Mass)
m, n = oh.mkline(mmin, 2.0, 1.0, 36.0)

# #############################################################################
# Loop/Plot
# #############################################################################
fig, ax = plt.subplots(1,1)
print "// Looping Outputs"
for nstep in nsteps:
    # Load
    print "// %012d/%012d" % (nstep, nsteps[-1])
    df = ioh.read_output("%s/Out_%s_%012d.dat" % (args.dir_name, \
                                                  args.run_name, \
                                                  nstep), \
                         frame="heliocentric")
    df = df[df.mass<16.0]
    s = df.mass * m + n
    dfx = df.sort(columns=["mass"], ascending=False).head(3)

    # Plot/Style
    ax.scatter(df.a, df.e, s=s**2.0, c="b", edgecolor="none", alpha=0.5)
    ax.set_title("%s/%s / %012d / %.2e yr" % (args.dir_name, \
                                              args.run_name, \
                                              nstep, \
                                              df.time.iloc[0]), \
                 fontsize="x-small")
    mtxt = "Most Massive Planets: (%.2f, %.2f, %.2f) M_Earth" % \
        ( dfx.iloc[0].mass, dfx.iloc[1].mass, dfx.iloc[2].mass )
    ax.text(0.05, 0.95, \
            mtxt, \
            horizontalalignment='left', \
            verticalalignment='top', \
            transform=ax.transAxes)
    ax.set_xlim([0,5])
    ax.set_ylim([0,0.6])
    ax.set_xlabel("Semi-Major Axis (AU)")
    ax.set_ylabel("Eccentricity")

    # Save
    fig.savefig("ae_%012d.png" % nstep)

    # Clean
    del df
    plt.cla()

del fig
print "// Done"
