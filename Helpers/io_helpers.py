"""
Pandas Helpers. Mostly to read Genga files.
"""

import numpy as np
import pandas as pd
import kepler_helpers as kh
import constants as C

# Single Genga Output
def read_output(fname):
    names_cols = [ "time", "pid", "mass", "radius", \
                   "x", "y", "z", \
                   "vx", "vy", "vz", \
                   "Sx", "Sy", "Sz", \
                   "amin", "amax", "emin", "emax", \
                   "aecount", "aecountT", "enccount", \
                   "test", "X" ]
    touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9 ]
    types_cols = { "pid": np.int32 }

    # Load CSV
    df = pd.read_csv(fname, \
                     sep=" ", \
                     header=None, \
                     names=names_cols, dtype=types_cols, \
                     usecols=touse_cols, \
                     index_col=1)

    # Convenience
    x = np.asarray(df.x); y = np.asarray(df.y); z = np.asarray(df.z)
    vx = np.asarray(df.vx); vy = np.asarray(df.vy); vz = np.asarray(df.vz)
    m = np.asarray(df.mass)

    # Barycentric Coordinates
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
def read_output_and_stack(fnames):
    names_cols = [ "time", "pid", "mass", "radius", \
                   "x", "y", "z", \
                   "vx", "vy", "vz", \
                   "Sx", "Sy", "Sz", \
                   "amin", "amax", "emin", "emax", \
                   "aecount", "aecountT", "enccount", \
                   "test", "X" ]
    touse_cols = [ 0, 1, 2, 4, 5, 6, 7, 8, 9 ]
    types_cols = { "pid": np.int32 }

    # Load CSV
    df = pd.DataFrame()
    for fname in fnames:
        try:
            dfx = pd.read_csv(fname, \
                              sep=" ", \
                              header=None, \
                              names=names_cols, dtype=types_cols, \
                              usecols=touse_cols, \
                              index_col=1)
            df = df.append(dfx)
        except IOError:
            pass

    # Drop Duplicate Indices
    # http://stackoverflow.com/questions/13035764/remove-rows-with-duplicate-indices-pandas-dataframe-and-timeseries
    df["index"] = df.index
    df.drop_duplicates(subset=["index"], inplace=True)
    del df["index"]

    # Convenience
    x = np.asarray(df.x); y = np.asarray(df.y); z = np.asarray(df.z)
    vx = np.asarray(df.vx); vy = np.asarray(df.vy); vz = np.asarray(df.vz)
    m = np.asarray(df.mass)

    # Barycentric Coordinates
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
def read_collisions_and_stack(fnames):
    names_cols = [ "time", \
                   "indexi", "mi", "ri", \
                   "xi", "yi", "zi", \
                   "vxi", "vyi", "vzi", \
                   "Sxi", "Syi", "Szi", \
                   "indexj", "mj", "rj", \
                   "xj", "yj", "zj", \
                   "vxj", "vyj", "vzj", \
                   "Sxj", "Syj", "Szj", "X" ]
    # touse_cols = [ 0, 1, 2, 4, 5, 6, 13, 14, 16, 17, 18 ]
    touse_cols = [ 0, 1, 2, 13, 14 ]
    types_cols = { "indexi": np.int32, "indexj": np.int32 }

    # Load CSV
    df = pd.DataFrame()
    for fname in fnames:
        dfx = pd.read_csv(fname, \
                          sep=" ", \
                          header=None, names=names_cols, dtype=types_cols, \
                          usecols=touse_cols)
        df = df.append(dfx, ignore_index=True)

    # Fix Mass
    dfx.mi *= C.msun/C.mearth
    dfx.mj *= C.msun/C.mearth

    # Return
    return df

