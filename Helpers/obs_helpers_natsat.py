"""
Observation Helpers for Natural Satellites.
http://localhost:9999/notebooks/_Thesis/02_SolarSystem/_Dev_JPL_HORIZONS_Moons.ipynb
"""

import numpy as np
import bs4 as bs4
import constants as C


def load_obs(fname='moons.html'):
    """
    Parse HTML Tables for Natural Satellites.
    Source: http://ssd.jpl.nasa.gov/?sat_phys_par

    @param: fname - HTML File w/ Moon Tables                [String]
    @return: mass - Moon Masses              (Earth Masses) [Numpy Float Array]
    """
    
    # Open/Parse HTML File.
    # May Fail if Format Changes!
    soup = bs4.BeautifulSoup(open(fname))
    tables = soup.find_all("table", attrs={"border": "1", \
                                           "cellspacing": "0", \
                                           "cellpadding": "5"})

    # Extract GM Columns
    GM = []
    for table in tables:
        rows = table.find_all("tr")[1:]
        for row in rows:
            GM_loc = float(row.find_all("td")[1].get_text().split(u'\xb1')[0])
            GM.append(GM_loc)

    # Extract Mass
    GM = np.asarray(GM) # km3/s2
    mass = GM/C.G       # kg
    mass /= C.mearth    # Mearth

    # Return
    return mass
