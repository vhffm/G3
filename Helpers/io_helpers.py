"""
I/O Helpers. Mostly to read Genga files.
"""

import numpy as np
import pandas as pd
import kepler_helpers as kh
import constants as C
import vector_helpers as vh

# Single Genga Output
def read_output(fname, frame):
    names_cols = [ "time", "pid", "mass", "radius", \
                   "x", "y", "z", \
                   "vx", "vy", "vz", \
                   "Sx", "Sy", "Sz", \
                   "amin", "amax", "emin", "emax", \
                   "aecount", "aecountT", "enccount", \
                   "test", "X" ]
    touse_cols = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 19 ]
    types_cols = { "pid": np.int32 }

    # User must pick reference frame
    # Genga outputs are heliocentric by default
    if not frame in [ "barycentric", "heliocentric" ]:
        estring = "Must Select Heliocentric/Barycentric Frame"
        raise Exception(estring)

    # Load CSV
    try:
        df = pd.read_csv(fname, \
                         sep=" ", \
                         header=None, \
                         names=names_cols, dtype=types_cols, \
                         usecols=touse_cols)
    except IOError:
         raise Exception("File Not Found: %s" % fname)

    # Convenience
    x = np.asarray(df.x); y = np.asarray(df.y); z = np.asarray(df.z)
    vx = np.asarray(df.vx); vy = np.asarray(df.vy); vz = np.asarray(df.vz)
    m = np.asarray(df.mass)

    # Barycentric Coordinates
    if frame == "barycentric":
        x, vx = kh.helio2bary(x, vx, m)
        y, vy = kh.helio2bary(y, vy, m)
        z, vz = kh.helio2bary(z, vz, m)

    # Compute Orbital Elements
    df["a"], df["e"], df["i"], df["Omega"], df["omega"], df["M"] = \
        kh.cart2kepX(x, y, z, vx, vy, vz, m)
    df["Q"] = df["a"] * ( 1 + df["e"] )
    df["q"] = df["a"] * ( 1 - df["e"] )

    # Fix Mass
    df.mass *= C.msun/C.mearth

    # Return
    return df

# Stack Multiple Genga Outputs, Remove Duplicate IDs
# fnames = [ fname01, fname02, ... ]
def read_output_and_stack(fnames, frame, drop_duplicates=True, nofail=False):
    names_cols = [ "time", "pid", "mass", "radius", \
                   "x", "y", "z", \
                   "vx", "vy", "vz", \
                   "Sx", "Sy", "Sz", \
                   "amin", "amax", "emin", "emax", \
                   "aecount", "aecountT", "enccount", \
                   "test", "X" ]
    touse_cols = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 19 ]
    types_cols = { "pid": np.int32 }

    # User must pick reference frame
    # Genga outputs are heliocentric by default
    if not frame in [ "barycentric", "heliocentric" ]:
        estring = "Must Select Heliocentric/Barycentric Frame"
        raise Exception(estring)

    # Load CSV
    df = pd.DataFrame()
    for ifname, fname in enumerate(fnames):
        try:
            dfx = pd.read_csv(fname, \
                              sep=" ", \
                              header=None, \
                              names=names_cols, dtype=types_cols, \
                              usecols=touse_cols)
            dfx['ifname'] = \
                pd.DataFrame({'ifname': np.ones(len(dfx)) * ifname})
            dfx['nstep'] = int(fname.split('/')[-1][:-4].split('_')[-1])
            df = df.append(dfx)
        except IOError:
            if not nofail:
                raise Exception("File Not Found: %s" % fname)

    # Drop Duplicate Indices
    # http://stackoverflow.com/questions/13035764/remove-rows-with-duplicate-indices-pandas-dataframe-and-timeseries
    if drop_duplicates:
        df.drop_duplicates(subset=["pid", "time", "nstep"], inplace=True)

    # Reindex (Relevant if we load multiple snapshots into one file)
    df.reset_index(drop=True, inplace=True)

    # Convenience
    x = np.asarray(df.x); y = np.asarray(df.y); z = np.asarray(df.z)
    vx = np.asarray(df.vx); vy = np.asarray(df.vy); vz = np.asarray(df.vz)
    m = np.asarray(df.mass)

    # Barycentric Coordinates
    if frame == "barycentric":
        x, vx = kh.helio2bary(x, vx, m)
        y, vy = kh.helio2bary(y, vy, m)
        z, vz = kh.helio2bary(z, vz, m)

    # Compute Orbital Elements
    df["a"], df["e"], df["i"], df["Omega"], df["omega"], df["M"] = \
        kh.cart2kepX(x, y, z, vx, vy, vz, m)
    df["Q"] = df["a"] * ( 1 + df["e"] )
    df["q"] = df["a"] * ( 1 - df["e"] )

    # Fix Mass
    df.mass *= C.msun/C.mearth

    # Return
    return df

# Stack Collision Files For Multiple Genga Outputs
# fnames = [ fname01, fname02, ... ]
def read_collisions_and_stack(fnames, return_xyz=False):
    names_cols = [ "time", \
                   "pidi", "mi", "ri", \
                   "xi", "yi", "zi", \
                   "vxi", "vyi", "vzi", \
                   "Sxi", "Syi", "Szi", \
                   "pidj", "mj", "rj", \
                   "xj", "yj", "zj", \
                   "vxj", "vyj", "vzj", \
                   "Sxj", "Syj", "Szj", "X" ]
    if return_xyz:
        touse_cols = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 19, 20, 21 ]
    else:
        touse_cols = [ 0, 1, 2, 13, 14 ]
    types_cols = { "pidi": np.int32, "pidj": np.int32 }

    # Load CSV
    df = pd.DataFrame()
    for ifname, fname in enumerate(fnames):
        try:
            dfx = pd.read_csv(fname, \
                              sep=" ", \
                              header=None, names=names_cols, \
                              dtype=types_cols, \
                              usecols=touse_cols)
            dfx['ifname'] = \
                pd.DataFrame({'ifname': np.ones(len(dfx)) * ifname})
            df = df.append(dfx, ignore_index=True)
        except IOError:
            raise Exception("File Not Found: %s" % fname)

    # Fix Mass
    df.mi *= C.msun/C.mearth
    df.mj *= C.msun/C.mearth

    # Fix Velocities, Return Relative Velocity (= Impact Velocity),
    #                        Angle between Velocities
    if return_xyz:
        df.vxi *= C.genga_to_kms
        df.vyi *= C.genga_to_kms
        df.vzi *= C.genga_to_kms
        df.vxj *= C.genga_to_kms
        df.vyj *= C.genga_to_kms
        df.vzj *= C.genga_to_kms
        df["dv"] = np.sqrt((df.vxi-df.vxj)**2.0 + \
                           (df.vyi-df.vyj)**2.0 + \
                           (df.vzi-df.vzj)**2.0)
        df["phi"] = vh.compute_angle(df.vxi, df.vyi, df.vzi, \
                                     df.vxj, df.vyj, df.vzj)

    # Return
    return df

# Stack Ejection Files For Multiple Genga Outputs
# fnames = [ fname01, fname02, ... ]
def read_ejections_and_stack(fnames):
    names_cols = [ "time", \
                   "pid", "m", "r", \
                   "x", "y", "z", \
                   "vx", "vy", "vz", \
                   "Sx", "Sy", "Sz", \
                   "case" ]
    # touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9, 13 ]
    touse_cols = [ 0, 1, 2, 13 ]
    types_cols = { "pid": np.int32, "case": np.int32 }

    # Load CSV
    df = pd.DataFrame()
    for ifname, fname in enumerate(fnames):
        try:
            dfx = pd.read_csv(fname, \
                              sep=" ", \
                              header=None, names=names_cols, \
                              dtype=types_cols, \
                              usecols=touse_cols)
            dfx['ifname'] = \
                pd.DataFrame({'ifname': np.ones(len(dfx)) * ifname})
            if len(dfx) > 0:
                df = df.append(dfx, ignore_index=True)
        except IOError:
            raise Exception("File Not Found: %s" % fname)

    # Empty?
    if len(df) == 0:
        df = pd.DataFrame({'time': [], \
                           'pid': [], 'm': [], 'case': [], \
                           'ifname': []})

    # Fix Mass
    df.m *= C.msun/C.mearth

    # Return
    return df


# Read (a,e,i) Grids
# Cf. http://localhost:9999/notebooks/HitnRun/Grid_AEI_Dev.ipynb
#
# Usage tips:
#
# ae_all[ae_all == 0] = np.nan
#
# fig, ax = plt.subplots(1,1)
# ax.imshow(ae_all, \
#           interpolation="None", cmap="Blues", origin="lower", \
#           extent=(0, 10, 0, 1))
# ax.set_aspect("auto")
def read_aei_grid(fname):

    # Load Data
    with open("%s" % fname, "r") as f:
        lines = f.readlines()

    #
    # Determine Blocks
    # Here, iendblock is the line after the last of the current block
    #

    iendblock = []
    first_emptyline = True
    for iline, line in enumerate(lines):
        if line.strip() == "":
            if first_emptyline:
                iendblock.append(iline)
                first_emptyline = False
            else:
                first_emptyline = True
    iendblock = np.asarray(iendblock, dtype=np.int64)

    #
    # Extract Matrices
    #

    # Reconstruct Dimensions
    Na = len(lines[iendblock[0]-1].strip().split(" "))
    Ne = iendblock[0]
    Ni = iendblock[2]-iendblock[1]-2

    # (a,e) 
    # All Time
    ae_all = np.ones([Ne, Na], dtype=np.int64)*np.nan
    for iline, line in enumerate(lines[:iendblock[0]]):
        ae_all[iline,:] = np.fromstring(line.strip(), sep=" ")

    # (a,e) 
    # Since Last Output
    ae_last = np.ones([Ne, Na], dtype=np.int64)*np.nan
    for iline, line in enumerate(lines[iendblock[0]+2:iendblock[1]]):
        ae_last[iline,:] = np.fromstring(line.strip(), sep=" ")
      
    # (a,i)
    # All Time
    ai_all = np.ones([Ni, Na], dtype=np.int64)*np.nan
    for iline, line in enumerate(lines[iendblock[1]+2:iendblock[2]]):
        ai_all[iline,:] = np.fromstring(line.strip(), sep=" ")

    # (a,i) 
    # Since Last Output
    ai_last = np.ones([Ni, Na], dtype=np.int64)*np.nan
    for iline, line in enumerate(lines[iendblock[2]+2:]):
        ai_last[iline,:] = np.fromstring(line.strip(), sep=" ")

    # Return
    return ae_all, ae_last, ai_all, ai_last
