"""
Remove moons from stdin (requires default Genga column ordering).

Call Signature:
$ python modic_cutmoons.py \
    --minpid 5 \
    --ce 95000 < $output_to_change > $changed_output
"""

import numpy as np
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--ce', type=int, required=True, \
                    help="Number of close encounters required for removal")
parser.add_argument('--minpid', type=int, default=1, \
                    help="Minimum PID that is cut.")
args = parser.parse_args()

# Read lines from stdin
lines_in = sys.stdin.read().rstrip("\n").split("\n")

# Remove particles with too many CEs
lines_out = []; cce = 0
for line in lines_in:
    pid = int(line.strip().split()[1])
    cnt = int(line.strip().split()[19])
    if (cnt < args.ce) or (pid < args.minpid):
        lines_out.append(line)
    else:
        sys.stderr.write("Removed Particle %i (%i CEs) \n" % \
            ( int(line.strip().split()[1]), int(line.strip().split()[19]) ) )
        cce += 1

# Removed Stats
sys.stderr.write("Removed %i Particles\n" % cce)

# Dump to stdout
for line in lines_out:
    print line
