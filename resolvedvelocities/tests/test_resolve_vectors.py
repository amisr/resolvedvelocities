# test_resolve_vectors.py

import sys
import os
sys.path.append(os.path.abspath('../'))
import ResolveVectors as rv


def main():
    config_file = sys.argv[1]
    vvels = rv.ResolveVectors(config_file)
    vvels.read_data()
    vvels.filter_data()
    vvels.transform()
    vvels.ion_upflow_correction()
    vvels.bin_data()
    vvels.get_integration_periods()
    vvels.compute_vector_velocity()
    vvels.compute_electric_field()
    vvels.compute_geodetic_output()
    vvels.save_output()

if __name__ == '__main__':
    main()
