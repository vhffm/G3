"""
Load Astorb File into Pandas Dataframe. Store as HDF5.
"""

import obs_helpers_minor as om
import pandas as pd
import argparse


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fname_in', default='astorb.dat', \
                    help='Name of Astorb Source File.')
parser.add_argument('-fname_out', default='astorb.hdf5', \
                    help='Name of Target HDF5 File.')
parser.add_argument('--short', action='store_true')
args = parser.parse_args()

# Info
if args.short:
    print "!! Short Mode. Loading 1000 Objects."
else:
    print "!! Loading 670k Objects."
    print "!! Get Some Coffee (Processing Time ~8 Minutes)."

# Process
print "// Loading Data (%s)" % args.fname_in
df = om.load_astorb(fname=args.fname_in, short=args.short)

# Store
print "// Saving to HDF5 (%s)" % args.fname_out
with pd.HDFStore(args.fname_out) as store:
    store["df"] = df
