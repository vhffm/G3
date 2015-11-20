"""
Extracts initial PIDs for sources.

Usage:
python ./extract_initial_pids.py --run_name run_01 > pid_list.csv

Returns one line of source PIDs per line in output. The first PID in the
list matches the PID of the final particle.
"""

import io_helpers as ioh
import formation_helpers as fh
import argparse
import sys
import numpy as np


# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-rname', '--run_name', \
                    default='run_01', \
                    help='Name of the Run.')
args = parser.parse_args()

# Load Collisions
fname = "Collisions_%s.dat" % args.run_name
dfc = ioh.read_collisions_and_stack([fname])

# Load IC Output
fname = "Out_%s_%012d.dat" % (args.run_name, 0)
dfo_t0 = ioh.read_output_and_stack([fname], frame='heliocentric')
dfo_t0 = dfo_t0[dfo_t0.mass < 12.0]

# Load Final Output
fname = "Out_%s_%012d.dat" % (args.run_name, 9e9)
dfo_tf = ioh.read_output_and_stack([fname], frame='heliocentric')
dfo_tf = dfo_tf[dfo_tf.mass < 12.0]

# Extract Final PIDs
final_pids = np.asarray(dfo_tf.pid, dtype=np.int64)
sys.stderr.write("// %i Final Particles\n" % len(final_pids))

# Loop Final PIDs
sources_all = []
for ipid, pid in enumerate(final_pids):
    sys.stderr.write("// Finding Sources for PID %i (%i/%i)\n" % \
        (pid, ipid+1, len(final_pids)))
    sources = fh.return_sources(pid, dfc)
    sources_all.append(sources)

# Output
for sources in sources_all:
    print ",".join(map(str, sources))
