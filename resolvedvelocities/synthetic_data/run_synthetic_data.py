import synthetic_data as synth
import numpy as np
import h5py

config_file_help = """Calculate 2D resolved plasma drift velocity and electric
field vectors from the LoS measurments in a fitted AMISR file."""

def main():

    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('synth_config_file',help='Configuration file for synthetic data set.')
    arg = parser.add_argument('vvels_config_file',help='Vvels config file.')
    args = vars(parser.parse_args())


    # synth_grid = np.meshgrid(np.linspace(62.,70.,10), np.linspace(260.,275.,10))
    # # synth_grid = np.meshgrid(np.linspace(82.,88.,10), np.linspace(305.,355.,10))
    # velocity = np.tile(np.array([500.,0.,0.]), (10,10,1))
    # synth_field = synth.Field(synth_grid[0], synth_grid[1], velocity, 300.)

    synth_field = synth.Field(args['synth_config_file'])


    # with h5py.File('/Volumes/AMISR_PROCESSED/processed_data/RISR-N/2017/11/imaginglp.v04.60.risrn/20171119.001/20171119.001_lp_1min-fitcal.h5', 'r') as f:
    #     beams = f['BeamCodes'][:,0]

    # az = [14.04,-154.30,-34.69,75.03]
    # el = [90.0, 77.5, 66.09, 65.56]
    # site = [65.13, -147.47, 0.213]

    # beams = np.array([64016, 64157, 64964, 65066])

    # radar = synth.Radar('PFISR', beams=beams)
    # # radar = synth.Radar(site, azimuth=az, elevation=el)
    radar = synth.Radar(args['synth_config_file'])

    synth.SyntheticData(synth_field, radar, args['vvels_config_file'])

if __name__ == '__main__':
    main()
