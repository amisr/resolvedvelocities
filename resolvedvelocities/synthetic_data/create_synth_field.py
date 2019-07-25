# create_synth_field.py
# script to print properly formated coordinate and field arrays for synth config file
# THIS SCRIPT DOES NOT GENERATE A FULL CONFIG FILE
# Instead, adust the parameters to get a desired field and run.  The output printed to the
#	screen can then be copy/pasted into a config file to avoid having to type out full arrays
# The synth config files were intentionally designed to just take in arrays of coordinates
#	and the field values at each point to be as flexible as possible in case future use wants
#	to evaluate multiple shears, vorticies, ect, but manually typing coordinates is a pain,
#	hence this script to help generate simple fields.

import numpy as np
from apexpy import Apex

glat = np.arange(62., 71., 1.)
glon = np.arange(200., 225., 3.)

glat, glon = np.meshgrid(glat, glon)
galt = np.full(glat.shape, 300.)
glat = glat.flatten()
glon = glon.flatten()
galt = galt.flatten()
coords = [glat.tolist(), glon.tolist(), galt.tolist()]

print('field_coords = ',coords)

# # Uniform field with only a Ve1 component
# Ve1 = 500.
# A = Apex(2019)
# f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(glat, glon, galt)
# Vgd = Ve1*e1.T
# Vgd = Vgd.tolist()

# Uniform geodetic East
Vgd = np.tile([500., 0., 0.], (len(galt),1)).tolist()

print('field_values = ', Vgd)