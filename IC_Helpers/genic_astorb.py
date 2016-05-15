"""
Generate Genga Initial Conditions File from Astorb Entries.

Usage:
$ python ./genic_astorb.py --fname_in astorb.hdf5 --crossers > initial.dat

Changes:

* 15 May 2016 / Volker Hoffmann <volker@cheleb.net>
  Add Perturbed Clone Generator (--clones)

* 15 May 2016 / Volker Hoffmann <volker@cheleb.net>
  Initial Release

"""

import sys
import pandas as pd
import kepler_helpers as kh
import constants as C
import argparse
import numpy as np


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fname_in', required=True, \
                    help='Name of Astorb Source File.')
parser.add_argument('--clones', action='store_true', \
                   help="Generate Perturbed Clones for Each Crosser.")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--crossers', action='store_true', \
                   help="Use Earth Crossers (~1'000).")
group.add_argument('--all', action='store_true', \
                   help="Use All Asteroid (~600'000).")
args = parser.parse_args()

# Load Data
sys.stderr.write('// Loading Data\n')
df = pd.read_hdf("%s" % args.fname_in, 'df')

# Select Earth Crossers
# Cf. ftp://cdsarc.u-strasbg.fr/pub/cats/B/astorb/ReadMe
if args.crossers:
    sys.stderr.write('// Select Earth Crossers\n')
    df = df[df.Xflg==1]

# Generate IC Lines
if args.clones:
    sys.stderr.write('// Generating Genga IC Lines, With Clones\n')
    pid_base = 10000; lines = []

    # Semi-Major Axis Offsets
    # a = a * [ 1,
    #           1-10^{-16}, 1-10^{-15}, ..., 1-10^{-1}
    #           1+10^{-1}, ..., 1+10^{-15}, 1-10^{-16} ]

    da_lo = 1.0 - np.logspace(-16, -1, 16)
    da_hi = 1.0 + np.logspace(-16, -1, 16)
    da_xx = np.concatenate([np.array([1.0]), da_lo, da_hi])

    for ii, [ irow, row ] in enumerate(df.iterrows()):
        for jj, da in enumerate(da_xx):
            # Convert to Cartesian State Vector
            x, v = \
                kh.kep2cart(row.a * da, row.e, row.i * C.d2r, \
                            row.Omega * C.d2r, row.omega * C.d2r, \
                            row.M * C.d2r, \
                            0.0, 1.0)
            v *= C.kms_to_genga

            # Format IC Line
            line = "0.0 %06d %.16e %.16e " % ( pid_base + jj, 0.0, \
                                               row.Diameter / C.au2km )
            line += "%+.16e %+.16e %+.16e " % ( x[0], x[1], x[2] )
            line += "%+.16e %+.16e %+.16e " % ( v[0], v[1], v[2] )
            line += "0.0 0.0 0.0"
            lines.append(line)

        # Base Particle ID Counter
        pid_base += 100

else:
    sys.stderr.write('// Generating Genga IC Lines\n')
    pid = 10000; lines = []
    for ii, [ irow, row ] in enumerate(df.iterrows()):
        # Convert to Cartesian State Vector
        x, v = \
            kh.kep2cart(row.a, row.e, row.i * C.d2r, \
                        row.Omega * C.d2r, row.omega * C.d2r, row.M * C.d2r, \
                        0.0, 1.0)
        v *= C.kms_to_genga

        # Format IC Line
        line = "0.0 %06d %.16e %.16e " % ( pid, 0.0, row.Diameter / C.au2km )
        line += "%+.16e %+.16e %+.16e " % ( x[0], x[1], x[2] )
        line += "%+.16e %+.16e %+.16e " % ( v[0], v[1], v[2] )
        line += "0.0 0.0 0.0"
        lines.append(line)

        # Particle IC Counter
        pid += 1

# Print Lines
sys.stderr.write('// Printing Genga IC Lines\n')
for line in lines:
    print line
