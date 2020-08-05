# run_resolve_vectors.py

import os
import sys
from .ResolveVectors import ResolveVectors

config_file_help = """Calculate 2D resolved plasma drift velocity and electric
field vectors from the LoS measurments in a fitted AMISR file.

Requires a configuration file containing the following example format:

[DEFAULT]

# Data file path/name
DATAFILE = /20190328.006/20190328.006_lp_1min-fitcal.h5

# output file name
OUTFILENAME = test_vvels.h5

# output path to save file in (optional)
OUTFILEPATH = 

# chirp
CHIRP = 0.0

# density limits - measurements with density outside of this range are discarded
NELIM = [2.0e9, 1.0e13]

# A priori covariance matrix
COVAR = [3000.*3000., 3000.*3000., 50.*50.]

# altitude limits - measurements outside of this range are discarded
ALTLIM = [150., 400.]

# chi2 limits - measurements with chi2 outside this range are discarded
CHI2LIM = [0.1, 10.]

# permitted fit codes - measurements with a fit code not in this list are discarded
GOODFITCODE = [1, 2, 3, 4]

# verticies of each bin in magnetic apex coordinates [alat, alon]
BINVERT = [[[65.0,-95.],[65.0,-88.],[65.5,-88.],[65.5,-95.]],
 [[65.5,-95.],[65.5,-88.],[66.0,-88.],[66.0,-95.]],
 [[66.0,-95.],[66.0,-88.],[66.5,-88.],[66.5,-95.]],
 [[66.5,-95.],[66.5,-88.],[67.0,-88.],[67.0,-95.]],
 [[67.0,-95.],[67.0,-88.],[67.5,-88.],[67.5,-95.]]]

# geodetic altitudes the output geodetic components should be calculated at
OUTALT = [100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.,850.]

# beamcode of the "up B" (oriented along the magnetic field) beam
UPB_BEAMCODE = 64157

# method for correcting for ion outflow
IONUP = UPB

# beamcodes of beams to be used (optional) - if omitted, all beams will be used
#USE_BEAMS = [64016, 64157, 64964]

# minumum number of points required to attemt to resolve a 2D vector
MINNUMPOINTS = 1

# post-integration time (optional) - if omitted, the native times of the input file are used
#INTTIME = 180.

# directory output summary plots should be saved in
PLOTSAVEDIR = /home/user/vvels/20190510.001_vvels_plots

"""


# a function to run this file as a script
def main():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('config_file',help='A configuration file.')

    args = vars(parser.parse_args())
    rv = ResolveVectors(args['config_file'])
    rv.read_data()
    rv.filter_data()
    rv.transform()
    rv.ion_upflow_correction()
    rv.bin_data()
    rv.get_integration_periods()
    rv.compute_vector_velocity()
    rv.compute_apex_velocity()
    rv.compute_electric_field()
    rv.compute_geodetic_output()
    rv.save_output()
    rv.create_plots()

if __name__=='__main__':
    main()
