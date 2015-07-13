"""
Extract (Stacked) Statistics (Min/Max/Avg/Std, 10/25/50/75/90 Quantiles):

- Mass
- Semi-Major Axis (Mass-Weighted)
- Eccentricity (Mass-Weighted)
- Inclination (Mass-Weighted)

for 

- All Particles
- Particles >= Cutoff
- Particles <  Cutoff

Usage: python /path/extract_stats2.py < dirlist

Here, "dirlist" is a file looking like:
/path/to/run_01
/path/to/run_02
...
/path/to/run_12
"""

import sys
import glob
import numpy as np
import pandas as pd
import io_helpers as ioh
import constants as C
import weighted as wq

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# CONFIG
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

m_cutoff = 2.0e23 # kg
m_cutoff /= C.mearth # earth masses

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# PREP
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

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

# The Globbit
run_names = []
for idir, cdir in enumerate(dirs):
    globs = glob.glob("%s/Out_*.dat" % cdir)
    # Extract run names
    # In:  Out_run_03_000057000000
    # Out: run_03
    run_names.append(globs[0].strip().split("/")[-1][:-4][4:-13])
    # Extract directory w/ most outputs
    if idir == 0:
        lmax = len(globs)
        xglobs = globs
    else:
        # Actual
        if len(globs) > lmax:
            lmax = len(globs)
            xglobs = globs

# Get steps
nsteps = np.zeros_like(xglobs, dtype=np.int64)
for iglob, xglob in enumerate(sorted(xglobs)):
    # In : /some/dir/Out_run_03_000156000000.dat
    # Out: 156000000
    nsteps[iglob] = int(xglob.strip().split("/")[-1][:-4].split("_")[-1])

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# ALLOCATE
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
# ALL PARTICLES
###############################################################################

# Time
time = np.zeros_like(nsteps, dtype=np.float64)

# Mass
mass_q10 = np.zeros_like(time)
mass_q25 = np.zeros_like(time)
mass_q50 = np.zeros_like(time)
mass_q75 = np.zeros_like(time)
mass_q90 = np.zeros_like(time)
mass_min = np.zeros_like(time)
mass_max = np.zeros_like(time)
mass_avg = np.zeros_like(time)
mass_std = np.zeros_like(time)

# Semi-Major Axis
a_q10 = np.zeros_like(time)
a_q25 = np.zeros_like(time)
a_q50 = np.zeros_like(time)
a_q75 = np.zeros_like(time)
a_q90 = np.zeros_like(time)
a_min = np.zeros_like(time)
a_max = np.zeros_like(time)
a_avg = np.zeros_like(time)
a_std = np.zeros_like(time)

# Eccentricity
e_q10 = np.zeros_like(time)
e_q25 = np.zeros_like(time)
e_q50 = np.zeros_like(time)
e_q75 = np.zeros_like(time)
e_q90 = np.zeros_like(time)
e_min = np.zeros_like(time)
e_max = np.zeros_like(time)
e_avg = np.zeros_like(time)
e_std = np.zeros_like(time)

# Inclination
i_q10 = np.zeros_like(time)
i_q25 = np.zeros_like(time)
i_q50 = np.zeros_like(time)
i_q75 = np.zeros_like(time)
i_q90 = np.zeros_like(time)
i_min = np.zeros_like(time)
i_max = np.zeros_like(time)
i_avg = np.zeros_like(time)
i_std = np.zeros_like(time)

###############################################################################
# ABOVE CUTOFF
###############################################################################

# Mass
mass_above_q10 = np.zeros_like(time)
mass_above_q25 = np.zeros_like(time)
mass_above_q50 = np.zeros_like(time)
mass_above_q75 = np.zeros_like(time)
mass_above_q90 = np.zeros_like(time)
mass_above_min = np.zeros_like(time)
mass_above_max = np.zeros_like(time)
mass_above_avg = np.zeros_like(time)
mass_above_std = np.zeros_like(time)

# Semi-Major Axis
a_above_q10 = np.zeros_like(time)
a_above_q25 = np.zeros_like(time)
a_above_q50 = np.zeros_like(time)
a_above_q75 = np.zeros_like(time)
a_above_q90 = np.zeros_like(time)
a_above_min = np.zeros_like(time)
a_above_max = np.zeros_like(time)
a_above_avg = np.zeros_like(time)
a_above_std = np.zeros_like(time)

# Eccentricity
e_above_q10 = np.zeros_like(time)
e_above_q25 = np.zeros_like(time)
e_above_q50 = np.zeros_like(time)
e_above_q75 = np.zeros_like(time)
e_above_q90 = np.zeros_like(time)
e_above_min = np.zeros_like(time)
e_above_max = np.zeros_like(time)
e_above_avg = np.zeros_like(time)
e_above_std = np.zeros_like(time)

# Inclination
i_above_q10 = np.zeros_like(time)
i_above_q25 = np.zeros_like(time)
i_above_q50 = np.zeros_like(time)
i_above_q75 = np.zeros_like(time)
i_above_q90 = np.zeros_like(time)
i_above_min = np.zeros_like(time)
i_above_max = np.zeros_like(time)
i_above_avg = np.zeros_like(time)
i_above_std = np.zeros_like(time)

###############################################################################
# BELOW CUTOFF
###############################################################################

# Mass
mass_below_q10 = np.zeros_like(time)
mass_below_q25 = np.zeros_like(time)
mass_below_q50 = np.zeros_like(time)
mass_below_q75 = np.zeros_like(time)
mass_below_q90 = np.zeros_like(time)
mass_below_min = np.zeros_like(time)
mass_below_max = np.zeros_like(time)
mass_below_avg = np.zeros_like(time)
mass_below_std = np.zeros_like(time)

# Semi-Major Axis
a_below_q10 = np.zeros_like(time)
a_below_q25 = np.zeros_like(time)
a_below_q50 = np.zeros_like(time)
a_below_q75 = np.zeros_like(time)
a_below_q90 = np.zeros_like(time)
a_below_min = np.zeros_like(time)
a_below_max = np.zeros_like(time)
a_below_avg = np.zeros_like(time)
a_below_std = np.zeros_like(time)

# Eccentricity
e_below_q10 = np.zeros_like(time)
e_below_q25 = np.zeros_like(time)
e_below_q50 = np.zeros_like(time)
e_below_q75 = np.zeros_like(time)
e_below_q90 = np.zeros_like(time)
e_below_min = np.zeros_like(time)
e_below_max = np.zeros_like(time)
e_below_avg = np.zeros_like(time)
e_below_std = np.zeros_like(time)

# Inclination
i_below_q10 = np.zeros_like(time)
i_below_q25 = np.zeros_like(time)
i_below_q50 = np.zeros_like(time)
i_below_q75 = np.zeros_like(time)
i_below_q90 = np.zeros_like(time)
i_below_min = np.zeros_like(time)
i_below_max = np.zeros_like(time)
i_below_avg = np.zeros_like(time)
i_below_std = np.zeros_like(time)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# LOOP STEPS
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

for istep, nstep in enumerate(nsteps):
    if istep % 32 == 0:
        print "// Processing %012d/%012d" % (nstep, nsteps[-1])
    fnames = []
    for idir, cdir in enumerate(dirs):
        fnames.append("%s/Out_%s_%012d.dat" % (cdir, run_names[idir], nstep))
    df = ioh.read_output_and_stack(fnames, \
                                   frame='heliocentric', \
                                   drop_duplicates=False, 
                                   nofail=True)
    df = df[df.mass < 12.0]
    df_above = df[df.mass >= m_cutoff]
    df_below = df[df.mass < m_cutoff]

    # Graceful Fail
    if len(df_above) == 0:
        sx = pd.Series(np.ones(len(df_above.columns.values)) * np.nan, \
                       index=df_above.columns.values)
        df_above = df_above.append(sx, ignore_index=True)

    # Time
    time[istep] = df.time.iloc[0]

    ###########################################################################
    # ALL PARTICLES
    ###########################################################################

    # Mass
    mass_q10[istep] = np.percentile(df.mass, 10)
    mass_q25[istep] = np.percentile(df.mass, 25)
    mass_q50[istep] = np.percentile(df.mass, 50)
    mass_q75[istep] = np.percentile(df.mass, 75)
    mass_q90[istep] = np.percentile(df.mass, 90)
    mass_min[istep] = df.mass.min()
    mass_max[istep] = df.mass.max()
    mass_avg[istep] = df.mass.mean()
    mass_std[istep] = df.mass.std()

    # Semi-Major Axis
    a_q10[istep] = wq.quantile_1D(df.a, df.mass, 0.10)
    a_q25[istep] = wq.quantile_1D(df.a, df.mass, 0.25)
    a_q50[istep] = wq.quantile_1D(df.a, df.mass, 0.50)
    a_q75[istep] = wq.quantile_1D(df.a, df.mass, 0.75)
    a_q90[istep] = wq.quantile_1D(df.a, df.mass, 0.90)
    a_min[istep] = df.a.min()
    a_max[istep] = df.a.max()
    a_avg[istep] = df.a.mean()
    a_std[istep] = df.a.std()

    # Eccentricity
    e_q10[istep] = wq.quantile_1D(df.e, df.mass, 0.10)
    e_q25[istep] = wq.quantile_1D(df.e, df.mass, 0.25)
    e_q50[istep] = wq.quantile_1D(df.e, df.mass, 0.50)
    e_q75[istep] = wq.quantile_1D(df.e, df.mass, 0.75)
    e_q90[istep] = wq.quantile_1D(df.e, df.mass, 0.90)
    e_min[istep] = df.e.min()
    e_max[istep] = df.e.max()
    e_avg[istep] = df.e.mean()
    e_std[istep] = df.e.std()

    # Inclination
    i_q10[istep] = wq.quantile_1D(df.i, df.mass, 0.10)
    i_q25[istep] = wq.quantile_1D(df.i, df.mass, 0.25)
    i_q50[istep] = wq.quantile_1D(df.i, df.mass, 0.50)
    i_q75[istep] = wq.quantile_1D(df.i, df.mass, 0.75)
    i_q90[istep] = wq.quantile_1D(df.i, df.mass, 0.90)
    i_min[istep] = df.i.min()
    i_max[istep] = df.i.max()
    i_avg[istep] = df.i.mean()
    i_std[istep] = df.i.std()

    ###########################################################################
    # ABOVE CUTOFF
    ###########################################################################

    # Mass
    mass_above_q10[istep] = np.percentile(df_above.mass, 10)
    mass_above_q25[istep] = np.percentile(df_above.mass, 25)
    mass_above_q50[istep] = np.percentile(df_above.mass, 50)
    mass_above_q75[istep] = np.percentile(df_above.mass, 75)
    mass_above_q90[istep] = np.percentile(df_above.mass, 90)
    mass_above_min[istep] = df_above.mass.min()
    mass_above_max[istep] = df_above.mass.max()
    mass_above_avg[istep] = df_above.mass.mean()
    mass_above_std[istep] = df_above.mass.std()

    # Semi-Major Axis
    a_above_q10[istep] = wq.quantile_1D(df_above.a, df_above.mass, 0.10)
    a_above_q25[istep] = wq.quantile_1D(df_above.a, df_above.mass, 0.25)
    a_above_q50[istep] = wq.quantile_1D(df_above.a, df_above.mass, 0.50)
    a_above_q75[istep] = wq.quantile_1D(df_above.a, df_above.mass, 0.75)
    a_above_q90[istep] = wq.quantile_1D(df_above.a, df_above.mass, 0.90)
    a_above_min[istep] = df_above.a.min()
    a_above_max[istep] = df_above.a.max()
    a_above_avg[istep] = df_above.a.mean()
    a_above_std[istep] = df_above.a.std()

    # Eccentricity
    e_above_q10[istep] = wq.quantile_1D(df_above.e, df_above.mass, 0.10)
    e_above_q25[istep] = wq.quantile_1D(df_above.e, df_above.mass, 0.25)
    e_above_q50[istep] = wq.quantile_1D(df_above.e, df_above.mass, 0.50)
    e_above_q75[istep] = wq.quantile_1D(df_above.e, df_above.mass, 0.75)
    e_above_q90[istep] = wq.quantile_1D(df_above.e, df_above.mass, 0.90)
    e_above_min[istep] = df_above.e.min()
    e_above_max[istep] = df_above.e.max()
    e_above_avg[istep] = df_above.e.mean()
    e_above_std[istep] = df_above.e.std()

    # Inclination
    i_above_q10[istep] = wq.quantile_1D(df_above.i, df_above.mass, 0.10)
    i_above_q25[istep] = wq.quantile_1D(df_above.i, df_above.mass, 0.25)
    i_above_q50[istep] = wq.quantile_1D(df_above.i, df_above.mass, 0.50)
    i_above_q75[istep] = wq.quantile_1D(df_above.i, df_above.mass, 0.75)
    i_above_q90[istep] = wq.quantile_1D(df_above.i, df_above.mass, 0.90)
    i_above_min[istep] = df_above.i.min()
    i_above_max[istep] = df_above.i.max()
    i_above_avg[istep] = df_above.i.mean()
    i_above_std[istep] = df_above.i.std()

    ###########################################################################
    # BELOW CUTOFF
    ###########################################################################

    # Mass
    mass_below_q10[istep] = np.percentile(df_below.mass, 10)
    mass_below_q25[istep] = np.percentile(df_below.mass, 25)
    mass_below_q50[istep] = np.percentile(df_below.mass, 50)
    mass_below_q75[istep] = np.percentile(df_below.mass, 75)
    mass_below_q90[istep] = np.percentile(df_below.mass, 90)
    mass_below_min[istep] = df_below.mass.min()
    mass_below_max[istep] = df_below.mass.max()
    mass_below_avg[istep] = df_below.mass.mean()
    mass_below_std[istep] = df_below.mass.std()

    # Semi-Major Axis
    a_below_q10[istep] = wq.quantile_1D(df_below.a, df_below.mass, 0.10)
    a_below_q25[istep] = wq.quantile_1D(df_below.a, df_below.mass, 0.25)
    a_below_q50[istep] = wq.quantile_1D(df_below.a, df_below.mass, 0.50)
    a_below_q75[istep] = wq.quantile_1D(df_below.a, df_below.mass, 0.75)
    a_below_q90[istep] = wq.quantile_1D(df_below.a, df_below.mass, 0.90) 
    a_below_min[istep] = df_below.a.min()
    a_below_max[istep] = df_below.a.max()
    a_below_avg[istep] = df_below.a.mean()
    a_below_std[istep] = df_below.a.std()

    # Eccentricity
    e_below_q10[istep] = wq.quantile_1D(df_below.e, df_below.mass, 0.10)
    e_below_q25[istep] = wq.quantile_1D(df_below.e, df_below.mass, 0.25)
    e_below_q50[istep] = wq.quantile_1D(df_below.e, df_below.mass, 0.50)
    e_below_q75[istep] = wq.quantile_1D(df_below.e, df_below.mass, 0.75)
    e_below_q90[istep] = wq.quantile_1D(df_below.e, df_below.mass, 0.90) 
    e_below_min[istep] = df_below.e.min()
    e_below_max[istep] = df_below.e.max()
    e_below_avg[istep] = df_below.e.mean()
    e_below_std[istep] = df_below.e.std()

    # Inclination
    i_below_q10[istep] = wq.quantile_1D(df_below.i, df_below.mass, 0.10)
    i_below_q25[istep] = wq.quantile_1D(df_below.i, df_below.mass, 0.25)
    i_below_q50[istep] = wq.quantile_1D(df_below.i, df_below.mass, 0.50)
    i_below_q75[istep] = wq.quantile_1D(df_below.i, df_below.mass, 0.75)
    i_below_q90[istep] = wq.quantile_1D(df_below.i, df_below.mass, 0.90) 
    i_below_min[istep] = df_below.i.min()
    i_below_max[istep] = df_below.i.max()
    i_below_avg[istep] = df_below.i.mean()
    i_below_std[istep] = df_below.i.std()

    del df, df_above, df_below

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# CREATE & STORE DATAFRAME
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Assemble dataframe panel
data = { \
    'time': time, \
    'mass_q10': mass_q10, \
    'mass_q25': mass_q25, \
    'mass_q50': mass_q50, \
    'mass_q75': mass_q75, \
    'mass_q90': mass_q90, \
    'mass_min': mass_min, \
    'mass_max': mass_max, \
    'mass_avg': mass_avg, \
    'mass_std': mass_std, \
    'a_q10': a_q10, \
    'a_q25': a_q25, \
    'a_q50': a_q50, \
    'a_q75': a_q75, \
    'a_q90': a_q90, \
    'a_min': a_min, \
    'a_max': a_max, \
    'a_avg': a_avg, \
    'a_std': a_std, \
    'e_q10': e_q10, \
    'e_q25': e_q25, \
    'e_q50': e_q50, \
    'e_q75': e_q75, \
    'e_q90': e_q90, \
    'e_min': e_min, \
    'e_max': e_max, \
    'e_avg': e_avg, \
    'e_std': e_std, \
    'i_q10': i_q10, \
    'i_q25': i_q25, \
    'i_q50': i_q50, \
    'i_q75': i_q75, \
    'i_q90': i_q90, \
    'i_min': i_min, \
    'i_max': i_max, \
    'i_avg': i_avg, \
    'i_std': i_std, \
    'mass_above_q10': mass_above_q10, \
    'mass_above_q25': mass_above_q25, \
    'mass_above_q50': mass_above_q50, \
    'mass_above_q75': mass_above_q75, \
    'mass_above_q90': mass_above_q90, \
    'mass_above_min': mass_above_min, \
    'mass_above_max': mass_above_max, \
    'mass_above_avg': mass_above_avg, \
    'mass_above_std': mass_above_std, \
    'a_above_q10': a_above_q10, \
    'a_above_q25': a_above_q25, \
    'a_above_q50': a_above_q50, \
    'a_above_q75': a_above_q75, \
    'a_above_q90': a_above_q90, \
    'a_above_min': a_above_min, \
    'a_above_max': a_above_max, \
    'a_above_avg': a_above_avg, \
    'a_above_std': a_above_std, \
    'e_above_q10': e_above_q10, \
    'e_above_q25': e_above_q25, \
    'e_above_q50': e_above_q50, \
    'e_above_q75': e_above_q75, \
    'e_above_q90': e_above_q90, \
    'e_above_min': e_above_min, \
    'e_above_max': e_above_max, \
    'e_above_avg': e_above_avg, \
    'e_above_std': e_above_std, \
    'i_above_q10': i_above_q10, \
    'i_above_q25': i_above_q25, \
    'i_above_q50': i_above_q50, \
    'i_above_q75': i_above_q75, \
    'i_above_q90': i_above_q90, \
    'i_above_min': i_above_min, \
    'i_above_max': i_above_max, \
    'i_above_avg': i_above_avg, \
    'i_above_std': i_above_std, \
    'mass_below_q10': mass_below_q10, \
    'mass_below_q25': mass_below_q25, \
    'mass_below_q50': mass_below_q50, \
    'mass_below_q75': mass_below_q75, \
    'mass_below_q90': mass_below_q90, \
    'mass_below_min': mass_below_min, \
    'mass_below_max': mass_below_max, \
    'mass_below_avg': mass_below_avg, \
    'mass_below_std': mass_below_std, \
    'a_below_q10': a_below_q10, \
    'a_below_q25': a_below_q25, \
    'a_below_q50': a_below_q50, \
    'a_below_q75': a_below_q75, \
    'a_below_q90': a_below_q90, \
    'a_below_min': a_below_min, \
    'a_below_max': a_below_max, \
    'a_below_avg': a_below_avg, \
    'a_below_std': a_below_std, \
    'e_below_q10': e_below_q10, \
    'e_below_q25': e_below_q25, \
    'e_below_q50': e_below_q50, \
    'e_below_q75': e_below_q75, \
    'e_below_q90': e_below_q90, \
    'e_below_min': e_below_min, \
    'e_below_max': e_below_max, \
    'e_below_avg': e_below_avg, \
    'e_below_std': e_below_std, \
    'i_below_q10': i_below_q10, \
    'i_below_q25': i_below_q25, \
    'i_below_q50': i_below_q50, \
    'i_below_q75': i_below_q75, \
    'i_below_q90': i_below_q90, \
    'i_below_min': i_below_min, \
    'i_below_max': i_below_max, \
    'i_below_avg': i_below_avg, \
    'i_below_std': i_below_std }

# Save
print "// Saving"
df = pd.DataFrame(data)
with pd.HDFStore("Stats2.hdf5", "w") as store:
    store['df'] = df
