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
    sources = np.asarray(sources, dtype=np.int64)
    return sources


def assign_wmf_raymond2004(a):
    """
    Cf. Sec. 2.2 and Fig. 2 in Raymond+ (2004).
    http://adsabs.harvard.edu/abs/2004Icar..168....1R

    @param: a - Semi Major Axis (Float)
    @return: wmf - Water Mass Fraction (Float)
    """
    # percent
    if a < 2.0:
        wmf = 0.001
    elif a >= 2.0 and a < 2.5:
        wmf = 0.1
    elif a>= 2.5:
        wmf = 5.0
    # fraction
    wmf /= 100.0
    # return
    return wmf


def assign_wmf_ronco2014(a):
    """
    Cf. Sec. 2 in Ronco+ (2014).
    http://adsabs.harvard.edu/abs/2015arXiv150907217R

    @param: a - Semi Major Axis (Float)
    @return: wmf - Water Mass Fraction (Float)
    """
    # percent
    if a < 2.0:
        wmf = 0.001
    elif a >= 2.0 and a < 2.5:
        wmf = 0.1
    elif a>= 2.5 and a < 2.7:
        wmf = 5.0
    else:
        wmf = 50.0
    # fraction
    wmf /= 100.0
    # return
    return wmf
