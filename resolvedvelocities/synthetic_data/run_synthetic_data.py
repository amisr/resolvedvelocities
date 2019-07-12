import synthetic_data as synth
import numpy as np

def main():
    synth_grid = np.meshgrid(np.linspace(62.,70.,10), np.linspace(260.,275.,10))
    velocity = np.tile(np.array([500.,0.,0.]), (10,10,1))
    synth_field = synth.Field(synth_grid[0], synth_grid[1], velocity, 300.)

    # az = [14.04,-154.30,-34.69,75.03]
    # el = [90.0, 77.5, 66.09, 65.56]
    site = [65.13, -147.47, 0.213]
    beams = [64016, 64157, 64964, 65066]
    radar = synth.Radar(site, beams=beams)
    # radar = synth.Radar(site, azimuth=az, elevation=el)

    synth.create_dataset(synth_field, radar)

if __name__ == '__main__':
    main()
