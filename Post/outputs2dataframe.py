"""
Convert Genga Coordinate Outputs to Dataframes (Per Directory).

Dirlist Format:
/path/01/
/path/02/
...
/path/NN/
"""

import io_helpers as ioh
import pandas as pd
import sys
import glob


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
for idir, cdir in enumerate(dirs):
    print "** %02d/%02d" % (idir, len(dirs))

    # Glob Coordinate Output Files
    globs = glob.glob("%s/Out_*.dat" % cdir)
    globs = sorted(globs)

    df = ioh.read_output_and_stack(globs, frame="heliocentric", \
                                   drop_duplicates=False)
    df.reset_index(drop=True, inplace=True)

    # Save
    with pd.HDFStore("Coordinates_%02d.hdf5" % int(idir+1)) as store:
        store["df"] = df

# Done
print "// Done"
