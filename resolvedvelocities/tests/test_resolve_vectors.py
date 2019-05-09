# test_resolve_vectors.py

import sys
import os
sys.path.append(os.path.abspath('../'))
import ResolveVectors as rv


def main():
    vvels = rv.ResolveVectors()
    vvels.read_data()
    vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()
    vvels.compute_vectors()

if __name__ == '__main__':
    main()
