import sys
import os
sys.path.append(os.path.abspath('../'))
# import synthetic_data as synth
from synthetic_data import Field, Radar
import numpy as np

def main():
    synth_grid = np.meshgrid(np.linspace(63.,69.,10), np.linspace(265.,270.,10))
    velocity = np.tile(np.array([300.,0.,0.]), (10,10,1))

    synth_field = Field(synth_grid[0], synth_grid[1], velocity, 300.)

    az = [14.04,-154.30,-34.69,75.03]
    el = [90.0, 77.5, 66.09, 65.56]
    site = [65.13, -147.47, 0.213]
    radar = Radar(site, az, el, 70.)

if __name__ == '__main__':
    main()
