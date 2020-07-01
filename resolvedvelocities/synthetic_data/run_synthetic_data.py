from resolvedvelocities.synthetic_data.SyntheticData import SyntheticData
from resolvedvelocities.synthetic_data.Field import Field
from resolvedvelocities.synthetic_data.Radar import Radar

config_file_help = """Create a synthetic AMISR dataset and calculate 2D resolved
plasma drift velocity vectors from the simulated LoS measurments.

Requires two configuration files:

Synthetic configuration file specifying the field and radar:

[[FIELD]

# year at which to initialize magnetic apex coordinate system
APEX_YEAR = 2019

# geodetic coordinates (glat, glon, galt) at which the field will be specified
FIELD_COORDS =  [[62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0], [200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 203.0, 203.0, 203.0, 203.0, 203.0, 203.0, 203.0, 203.0, 203.0, 206.0, 206.0, 206.0, 206.0, 206.0, 206.0, 206.0, 206.0, 206.0, 209.0, 209.0, 209.0, 209.0, 209.0, 209.0, 209.0, 209.0, 209.0, 212.0, 212.0, 212.0, 212.0, 212.0, 212.0, 212.0, 212.0, 212.0, 215.0, 215.0, 215.0, 215.0, 215.0, 215.0, 215.0, 215.0, 215.0, 218.0, 218.0, 218.0, 218.0, 218.0, 218.0, 218.0, 218.0, 218.0, 221.0, 221.0, 221.0, 221.0, 221.0, 221.0, 221.0, 221.0, 221.0, 224.0, 224.0, 224.0, 224.0, 224.0, 224.0, 224.0, 224.0, 224.0], [300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0]]

# geodetic components (east, north, up) of the velocity field at each coordinate
FIELD_VALUES =  [[500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0], [500.0, 0.0, 0.0]]


[RADAR]

# radar site coordinates (geodetic lat, lon, alt)
SITE_COORDS = [65.13,-147.47,0.213]

# file name for the beam code table
BEAMCODE_FILENAME = 'bcotable_pfisr.txt'

# beam codes to use
BEAMCODES = [64016, 64157, 64964, 65066]

# lenght of each range gate
RANGE_STEP = 50.

# name of output file containing the synthetic data set - this file format mimics the standard AMISR fitted file format and will be used in the normal resolvedvelocities algorithm
OUTPUT_FILENAME = 'synthetic_data.h5'

# integration period of each time step in simulated data set (this doesn't really matter)
INTEGRATION_PERIOD = 60.

# measurment error to assign to each LoS velocity
VEL_ERROR = 10.



Standard vvels configuration file containing the following example format:

[DEFAULT]

# Data file path/name
DATAFILE = synthetic_data.h5

# output file name
OUTFILENAME = test_vvels.h5

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
"""

def main():

    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('synth_config_file',help='Configuration file for synthetic data set.')
    arg = parser.add_argument('vvels_config_file',help='Vvels config file.')
    args = vars(parser.parse_args())


    # generate field object
    field = Field(args['synth_config_file'])
    # field.plot_ionosphere()

    # generate radar object
    radar = Radar(args['synth_config_file'])

    # use field and radar objects to produce synthetic data set
    synth_data = SyntheticData(field, radar)

    # run resolvevectors algothrithm on synthetic data set
    vvels = synth_data.eval_vvels(args['vvels_config_file'])

    # compare output of resolvevectors algorithm with truth
    synth_data.plot(field, radar, vvels)
    # synth_data.check_assumptions(field, radar, vvels)
    synth_data.check_components(field, vvels)


if __name__ == '__main__':
    main()
