# test_resolve_vectors.py

import sys
import os
sys.path.append(os.path.abspath('../'))
import ResolveVectors as rv


def main():
    vvels = rv.ResolveVectors()
    vvels.read_data()
    vvels.filter_data()

if __name__ == '__main__':
    main()
