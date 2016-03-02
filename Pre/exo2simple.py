"""
Convert Exoplanet.Eu Dump to Simple CSV Format.
For Gnuplot Lovers and Python Haters.
Get Dump Here: http://exoplanet.eu/catalog/csv/
"""

import pandas as pd

# Constants
mearth   = 5.97219e24 # kg
mjupiter = 1.89813e27 # kg
Rearth   = 6371.0  # km
Rjupiter = 69911.0 # km

# Source
# fexo = '/home/ics/volker/ExoEu/exoplanet.eu_catalog_2016-01-11.csv'
fexo = 'exoplanet.csv'

# Load
print "// Loading"
df = pd.read_csv(fexo, header=0)
df.rename(columns = {'# name': 'name'}, inplace = True)
# df.mass = df.mass * mjupiter / mearth
# df.radius = df.radius * Rjupiter / Rearth

# Print Some Data
print "// SAMPLE OUTPUT"
print df[['name','star_name','radius','mass',\
          'semi_major_axis','eccentricity']].head(3).T

# Export
# print "// Saving"
# df.to_csv('exoplanet_simple.csv', \
#           columns=['mass', 'semi_major_axis', 'radius'], \
#           index=False)

# Use this if you want ? in unnamed columns
print "// Saving"
df.to_csv('exoplanet_simple.csv', \
          columns=['mass', 'semi_major_axis', 'radius'], \
          na_rep='?', \
          index=False)