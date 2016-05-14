"""
Observation Helpers for Minor Bodies.
http://localhost:9999/notebooks/Astorb/Astorb.ipynb
"""

import pandas as pd
import numpy as np


def load_astorb(fname='astorb.dat', short=False):
    """
    Load Astorb Database from Fixed-Width File.
    Includes Data for Minor Bodies. Excludes Natural Satellites.
    
    Raw Data : Orbital Elements, Magnitudes, Types, Albedos, ...
    Derived  : Type => Albedo, Density
             : B-V Mag + Albedo => Diameter
             : Diameter + Density => Mass

    Conversions from Caption of Fig. 1 in Petit+ 2001
    Cf. http://adsabs.harvard.edu/abs/2001Icar..153..338P
    Cf. http://www.minorplanetcenter.net/iau/lists/Sizes.html

    Takes ~8 Minutes for Full File (~670k Rows).
    Cf. ftp://cdsarc.u-strasbg.fr/pub/cats/B/astorb/astorb.html

    Details on Columns in ReadMe.
    Cf. ftp://cdsarc.u-strasbg.fr/pub/cats/B/astorb/ReadMe

    @param: fname - Filename for Asteroid Database [String]
    @return: df   - Minor Bodies Orbital Databse   [Pandas Dataframe]
    """

    # ftp://cdsarc.u-strasbg.fr/pub/cats/B/astorb/astorb.html
    # http://pandas.pydata.org/pandas-docs/dev/io.html#files-with-fixed-width-columns
    # so, first number is one below/before the col. number of the first letter
    # end is the column number of the last letter

    # Ka-Ching
    colspecs = [ (0, 6), (7, 25), (26, 41), (42, 48), \
                 (49, 53), (54, 58), (59, 64), (65, 72), \
                 (71, 74), (75, 78), (79, 82), (83, 86), (87, 90), (91, 94), \
                 (95, 100), (101, 105), (106, 114), (115, 125), (126, 136), \
                 (137, 147), (148, 157), (158, 168), (169, 181), \
                 (182, 190), (191, 198), (199, 207), (208, 216), \
                 (217, 224), (225, 233), (234, 241), (242, 250), \
                 (251, 258), (259, 267) ]

    # H = Absolute Magnitude
    # B-V = Colour Index
    colnames = [ "Number", "Name", "Computer", "Magnitude", \
                 "Slope", "B-V", "Diameter", "Type", \
                 "Xflg", "Oflg", "Sflg", "MPCC", "LowC", "FlaC", \
                 "Arc", "NObservations", "Epoch", "M", "omega", "Omega", \
                 "i", "e", "a", "ComputeDate", "CEU", "CEU-Rate", "CEU-Date", \
                 "PEU", "PEU-Max", "PEU-Max10", "PEU-Max10-Date", \
                 "PEX-Max10-X", "PEU-Max10-X-Date" ]

    # DO NOT USE THE ID TO INDEX; IT MAY BE BLANK
    # df = pd.read_fwf(fname, colspecs=colspecs, header=None, index_col=0)

    # Short Mode for Debugging?
    if short:
        df = pd.read_fwf(fname, colspecs=colspecs, \
                         header=None, names=colnames, \
                         parse_dates=[16,23,26,28,30,32], nrows=1000)
    else:
        df = pd.read_fwf(fname, colspecs=colspecs, header=None, \
                         names=colnames, \
                         parse_dates=[16,23,26,28,30,32])

    # Compute qQ
    df["Q"] = df["a"] * ( 1.0 + df["e"] )
    df["q"] = df["a"] * ( 1.0 - df["e"] )

    # Fix Strings
    # http://pandas.pydata.org/pandas-docs/dev/text.html
    df["Type"] = df["Type"].str.upper()
    df["Type"] = df["Type"].str.replace('?', '')
    df["Type"] = df["Type"].str[0]

    #
    # THE FOLLOWING ONLY HOLDS FOR THE ASTEROID BELT
    #

    # Albedo Table
    types = [ "K", "M", "S", "A", "E", "C", "D", "F", \
              "G", "R", "V", "L", "P", "T", np.nan ]
    albedos = np.array([ 0.20, 0.20, 0.20, 0.40, 0.40, 0.05, 0.05, 0.05, \
                         0.05, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 ])
    data_albedo = { "Type": types, "Albedo": albedos } 
    df_albedo = pd.DataFrame(data=data_albedo)
    df_albedo[["Type", "Albedo"]].T

    # Density Table
    types = [ "C", "D", "F", "G", "K", "A", "E", "S", \
              "M", "R", "V", "L", "P", "T", np.nan ]
    densities = np.array([ 1.5, 1.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5, \
                           3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 ])
    data_density = { "Type": types, "Density": densities }
    df_density = pd.DataFrame(data=data_density)
    df_density[["Type", "Density"]].T

    # Merge Properties
    df_properties = df_albedo.merge(df_density, how="left", on="Type")
    # df_properties[["Type", "Albedo", "Density"]].T
    # df_properties.T

    # Merge Albedo, Density
    df = df.merge(df_properties, how="left", on="Type")

    #
    # ADJUST FOR OUTER SOLAR SYSTEM
    #
    # Known Albedos: 
    # - Eris (0.96)
    # - Pluto (0.49 to 0.66)
    # - Charon (0.2 to 0.6)
    # - 67p (0.06) [Thrown In Comet, Fluffy!]
    # - Icy Stuff (0.5; http://www.minorplanetcenter.net/iau/lists/Sizes.html)
    # - Chaos (0.04 to 0.08)
    # - Varuna (0.08 to 0.16)
    # - Sedna (0.26 to 0.38)
    #
    # Known Densities (g/cm3)
    # - Eris (2.52)
    # - Pluto (1.86)
    # - Charon (1.7)
    # - 67p (0.4) [Throw In Comet, Fluffy!]
    # - Varuna (0.992)
    #
    # For known albedos/densities that we bothered to tabulate, use them.
    # Otherwise, use 0.6 (albedo), and 1.2 g/cm3 (density)
    #

    # Catch-All
    df.loc[df.a > 4.5, 'Albedo' ] = 0.6
    df.loc[df.a > 4.5, 'Density'] = 1.2

    # Tabulated Objects
    # Albedo
    df.loc[df.Name=='Eris',     'Albedo']  = 0.96
    df.loc[df.Name=='Pluto',    'Albedo']  = 0.55
    df.loc[df.Name=='Charon',   'Albedo']  = 0.40
    df.loc[df.Name=='Chaos',    'Albedo']  = 0.06
    df.loc[df.Name=='Varuna',   'Albedo']  = 0.12
    df.loc[df.Name=='Sedna',    'Albedo']  = 0.32
    df.loc[df.Name=='Makemake', 'Albedo']  = 0.80
    df.loc[df.Name=='Haumea',   'Albedo']  = 0.75

    # Density
    df.loc[df.Name=='Eris',     'Density'] = 2.52
    df.loc[df.Name=='Pluto',    'Density'] = 1.86
    df.loc[df.Name=='Charon',   'Density'] = 1.70
    df.loc[df.Name=='Varuna',   'Density'] = 0.99
    df.loc[df.Name=='Makemake', 'Density'] = 2.30
    df.loc[df.Name=='Haumea',   'Density'] = 2.95

    # Compute Missing Diameters
    # @todo - Ref???
    idx = df.Diameter.isnull()
    df.loc[idx,"Diameter"] = \
        10.0**(0.5*(6.244 - \
                    0.4 * df.loc[idx,"Magnitude"] - \
                    np.log10(df.loc[idx,"Albedo"])))

    # Compute Masses
    df["Mass"] = (df["Diameter"]/2.0)**3.0 * 4./3. * np.pi * \
                 (df["Density"]*100.0**3.0*1000.0**3.0/1000.0)
    # df[["Name", "Type", "Magnitude", "Diameter", "Albedo", "Density", "Mass"]].head(10)

    # Return
    return df


def load_mpcorb(fname='PHA.txt'):
    """
    http://www.minorplanetcenter.net/iau/MPCORB.html
    http://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html
    """
