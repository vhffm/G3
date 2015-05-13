"""
Helpers for Formation Simulations.
"""

import numpy as np


def return_sources(pid, dfc):
    """
    Construct List of Mass Sources for a given Particle ID.
    Can be used to build Merger Trees.

    @param pid: Particle ID - [Integer]
    @param dfc: Collision list - [Pandas Dataframe from Io_Helpers] 
    @return sources: Source Particle IDs - [Numpy Array]
    """
    dfc = dfc.sort(columns="time", ascending=True)
    sources = [pid]
    for irow in range(len(dfc)):
        irow = len(dfc)-irow-1
        dfc_loc = dfc.iloc[irow]
        if dfc_loc.pidi in sources:
            sources.append(int(dfc_loc.pidj))
        elif dfc_loc.pidj in sources:
            sources.append(int(dfc_loc.pidi))
    sources = np.asarray(source, dtype=np.int64)
    return sources
