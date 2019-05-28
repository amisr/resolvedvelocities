import sys
import os
sys.path.append(os.path.abspath('../'))
# import synthetic_data as synth
from synthetic_data import Field
import numpy as np

def main():
    synth_grid = np.meshgrid(np.linspace(63.,69.,10), np.linspace(265.,270.,10))
    velocity = np.tile(np.array([300.,0.,0.]), (10,10,1))

    synth_field = Field(synth_grid[0], synth_grid[1], velocity, 300.)


if __name__ == '__main__':
    main()
