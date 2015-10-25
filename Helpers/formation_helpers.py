"""
Helpers for Formation Simulations.
"""

import numpy as np
import stats_helpers as sh
import pandas as pd


def return_sources(pid, dfc):
    """
    Construct List of Mass Sources for a given Particle ID.
    Can be used to build Merger Trees.

    @todo: Accelerate. Rewrite in Fortran?

    @param pid: Particle ID - [Integer]
    @param dfc: Collision list - [Pandas Dataframe from Io_Helpers] 
    @return sources: Source Particle IDs - [Numpy Array]
    """
    dfc = dfc.sort(columns="time", ascending=True)
    sources = np.zeros(len(dfc)+1) * np.nan
    sources[0] = pid
    for ii, irow in enumerate(range(len(dfc))):
        irow = len(dfc)-irow-1
        dfc_loc = dfc.iloc[irow]
        if dfc_loc.pidi in sources:
            sources[ii+1] = int(dfc_loc.pidj)
        elif dfc_loc.pidj in sources:
            sources[ii+1] = int(dfc_loc.pidi)
    sources = sources[~np.isnan(sources)]
    sources = np.asarray(source, dtype=np.int64)
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
    http://adsabs.harvard.edu/abs/2014A&A...567A..54R

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


def compute_kde(df, evaluation_range, evaluation_range_step, variable, \
                cov_tight_factor=4.0):
    """
    Compute Kernel Density Estimate for Formation Runs.

    Options:
    1. Mass-Weighted Semi-Major Axis Distribution
    2. Mass-Weighted Eccentricity Distribution
    3. Mass-Weighted Inclination Distribution
    4. Mass Function.

    Picking a cov_tight_factor is a bit of a dark art. Too small and the 
    KDE is way too smooth. Too large and the resulting KDE doesn't smooth over
    adjacent points at all. Usually, power of 2 in the range 2 to 12 are 
    reasonable. Use common sense. (@todo: Add literature link.)

    @param: df - Genga Coordinate DataFrame [Pandas DataFrame]
    @param: evaluation_range - Range over which to compute KDE [Np Float Array]
    @param: evaluation_range_step - Step size for Evaluation Range [Float]
    @param: variable: - Variable for which to compute KDE (a,e,i,m) [String]
    @param: cov_tight_factor - Tighten covariance by this factor [Float]
    @return: kde_evaluated - KDE evaluated on evaluation_range [Np Float Array]
    """

    # Select Variable to Fit KDE
    if variable == 'a':
        input_array = np.asarray(df.a)
    elif variable == 'e':
        input_array = np.asarray(df.e)
    elif variable == 'i':
        input_array = np.asarray(df.i)
    elif variable == 'm':
        input_array = np.log10(np.asarray(df.mass))
    else:
        raise Exception("Invalid Variable Requested.")

    # Set Weights
    if variable in [ 'a', 'e', 'i' ]:
        weights = np.asarray(df.mass/df.mass.max())
    else:
        weights = np.ones_like(df.mass)
        
    # Compute KDE
    cv = sh.Covariator(np.atleast_2d(input_array), weights)
    inv_cov, norm_factor = cv(bw_method = 'scott')
    inv_cov *= cov_tight_factor
    kde = sh.gaussian_kde(np.atleast_2d(input_array), weights, \
                       inv_cov, norm_factor)
    kde_evaluated = kde(evaluation_range)

    # Normalize to Total Mass (Int{dM/d{a,e,i}} = Mass)
    if variable in [ 'a', 'e', 'i' ]:
        norm = np.sum(df.mass) / ( np.sum(kde_evaluated) * evaluation_range_step )

    # Normalize to Total Number of Particles (Int{dN/dM} = N)
    else:
        dmx = np.diff(10.0**evaluation_range) # Midpoints in mrange, convert to lin
        kde_evaluated_x = \
            (kde_evaluated[1:] + kde_evaluated[:-1]) / 2.0 # Linear Interp.
        N_0 = np.sum(kde_evaluated_x * dmx)
        norm = len(df) / N_0

    # Apply Normalization
    kde_evaluated *= norm

    # Return
    return kde_evaluated


def compute_wmf(dfo, dfo_t0, dfc, showstep=False):
    """
    Compute Water Mass Fraction for all Particles.
    Build Source List for Output, Compute WMF from Precursors.

    @param: dfo - Coordinate Output @ Time [Pandas Dataframe]
    @param: dfo_t0 - Outout @ Initial Time [Pandas Dataframe]
    @param: dfc - Collision List [Pandas Dataframe]
    @return: dfo - Coordinate Output @ Time w/ WMF Fields [Pandas Dataframe]
    """

    wmf_01 = np.ones(len(dfo)) * np.nan
    wmf_02 = np.ones(len(dfo)) * np.nan
    
    dfc_now = dfc[dfc.time<=dfo.time.iloc[0]]
    
    # The Murder Loop
    # @todo - Accelerate? Rewrite Source List Construction in Fortran?
    for ii, dfo_row in enumerate(dfo.iterrows()):
        # Info
        if showstep:
            if ii % int(len(dfo)/8) == 0:
                print "%i/%i" % (ii,len(dfo))
        
        # Extract Series from One-Row Dataframe
        dfo_loc = dfo_row[1]
        
        # Identify Source Particles
        sources = return_sources(int(dfo_loc.pid), dfc_now)
        dfo_sources = dfo_t0[dfo_t0.pid.isin(sources)]
        
        # Compute WMF (Raymond+ 2004, Ronco+ 2014)
        dfo_sources.loc[:,'wmf_01'] = \
            pd.Series(dfo_sources.a.apply(assign_wmf_raymond2004), \
                      index=dfo_sources.index)
        dfo_sources.loc[:,'wmf_02'] = \
            pd.Series(dfo_sources.a.apply(assign_wmf_ronco2014), \
                      index=dfo_sources.index)
        wmf_01[ii] = np.sum(dfo_sources.mass * dfo_sources.wmf_01) / \
            np.sum(dfo_sources.mass)
        wmf_02[ii] = np.sum(dfo_sources.mass * dfo_sources.wmf_02) / \
            np.sum(dfo_sources.mass)
    
    # Append WMF
    dfo.loc[:,'wmf_01'] = wmf_01
    dfo.loc[:,'wmf_02'] = wmf_02
    
    # Return
    return dfo


def compute_wmf_wrapper(wrapper):
    """
    Wrapper for WMF Computation.
    Expects Arguments in List. Unfolds, and Calls WMF Computation.

    @param: wrapper - List of Coordinates, Collisions [Python List]
                      Cf. compute_wrapper() Doc
    @return: dfo - Coordinate Output @ Time w/ WMF Fields [Pandas Dataframe]
    """
    
    # Unwrap
    dfo = wrapper[0]
    dfo_t0 = wrapper[1]
    dfc = wrapper[2]
    
    # Call Unwrapped
    dfo = compute_wmf(dfo, dfo_t0, dfc)
    
    # Return
    return dfo
