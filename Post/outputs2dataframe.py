"""
Convert Genga Coordinate Outputs to Dataframes (Per Directory).

Dirlist Format:
/path/01/
/path/02/
...
/path/NN/

Code is embarrassingly parallel. Use the -np to define number of subprocesses.
"""

import io_helpers as ioh
import pandas as pd
import sys
import glob
import multiprocessing as mp
import argparse


###############################################################################
# FUNCTION DEFINITIONS
###############################################################################

def outs2df(cdir):
    """
    Processing Function.
    """

    print "** %s" % cdir

    # Glob Coordinate Output Files
    globs = glob.glob("%s/Out_*.dat" % cdir)
    globs = sorted(globs)

    # Load Data
    df = ioh.read_output_and_stack(globs, frame="heliocentric", \
                                   drop_duplicates=False)
    df.reset_index(drop=True, inplace=True)

    # Return Dataframe
    return df


###############################################################################
# MAIN PROGRAM STARTS HERE
###############################################################################

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-np', type=int, default=1, \
                    help='Number of Processes')
args = parser.parse_args()
print "// Using %i Subprocesses" % args.np

# List of Directories
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    for line in lines:
        dirs.append(line)
    print "// Reading %i Directories" % len(dirs)


# Loop Directories
if args.np == 1:
    result = []
    for idir, cdir in enumerate(dirs):
        df_tmp = outs2df(cdir)
        result.append(df_tmp)
else:
    pool = mp.Pool(processes=args.np)
    result = pool.map(outs2df, dirs)
    pool.close()
    pool.join()

# Save
print "// Saving Dataframes"
for idf, df in enumerate(result):
    with pd.HDFStore("Coordinates_%02d.hdf5" % int(idf+1)) as store:
        store["df"] = df

# Done
print "// Done"
