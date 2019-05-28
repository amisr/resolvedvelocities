import sys
import os
sys.path.append(os.path.abspath('../'))
import synthetic_data as synth
import numpy as np

def main():
    synth_grid = np.meshgrid(np.linspace(63.,69.,50), np.linspace(265.,270.,50))
    velocity = np.tile(np.array([300.,0.,0.]), (50,50,1))

    synth.map_velocity_field(synth_grid, velocity, 300.)

if __name__ == '__main__':
    main()
