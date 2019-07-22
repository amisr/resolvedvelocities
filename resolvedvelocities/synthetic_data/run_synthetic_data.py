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


    # generate field object
    synth_field = synth.Field(args['synth_config_file'])

    # generate radar object
    radar = synth.Radar(args['synth_config_file'])

    # use field and radar objects to produce synthetic data set
    synth_data = synth.SyntheticData(synth_field, radar)

    # run resolvevectors algothrithm on synthetic data set
    vvels = synth_data.eval_vvels(args['vvels_config_file'])

    # compare output of resolvevectors algorithm with truth
    synth_data.plot(synth_field, radar, vvels)
    synth_data.compare_components(synth_field,vvels)


if __name__ == '__main__':
    main()
