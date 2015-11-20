"""
Compute Collision Geometries and Energies.
Processing ~1900 Collisions Takes ~ 1min30s.
"""

import io_helpers as ioh
import formation_helpers as fh
import collision_helpers as ch
import pandas as pd
import argparse
import sys
import glob
import numpy as np


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fin', '--input_file', \
                    default='Collisions_run_01.dat', \
                    help='Name of Output File.')
parser.add_argument('-fout', '--output_file', default='Collisions.hdf5', \
                    help='Name of Output File.')
args = parser.parse_args()

# Load Collisions
print "// Loading Collisions"
dfc = ioh.read_collisions_and_stack([args.input_file], return_xyz=True)

# Reconstruct Geometries
print "// Reconstructing %i Geometries" % len(dfc)
theta, b_over_r, v_impact = ch.reconstruct_geometries(dfc.copy())
dfc['theta'] = theta
dfc['b_over_r'] = b_over_r
dfc['v_impact'] = v_impact
nfailed = np.sum(np.isnan(dfc.theta))
print "// Failed to Reconstruct %i/%i Geometries (%.2f%%)." % \
    ( nfailed , len(dfc), round((float(nfailed)/float(len(dfc)) * 100.0)) )

# Compute Specific Impact Energies
print "// Computing Impact Energies"
dfc = ch.compute_specific_impact_energy(dfc)

# Save
print "// Saving to %s" % args.output_file
with pd.HDFStore("%s" % args.output_file, 'w') as store:
    store['df'] = dfc
