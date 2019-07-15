import synthetic_data as synth
import numpy as np
import h5py

def main():
    synth_grid = np.meshgrid(np.linspace(62.,70.,10), np.linspace(260.,275.,10))
    velocity = np.tile(np.array([500.,0.,0.]), (10,10,1))
    synth_field = synth.Field(synth_grid[0], synth_grid[1], velocity, 300.)


    with h5py.File('/Volumes/AMISR_PROCESSED/processed_data/RISR-N/2017/11/imaginglp.v04.60.risrn/20171119.001/20171119.001_lp_1min-fitcal.h5', 'r') as f:
        beams = f['BeamCodes'][:,0]

    # az = [14.04,-154.30,-34.69,75.03]
    # el = [90.0, 77.5, 66.09, 65.56]
    # site = [65.13, -147.47, 0.213]

    # beams = np.array([64016, 64157, 64964, 65066])
    print(beams)

    radar = synth.Radar('RISRN', beams=beams)
    # radar = synth.Radar(site, azimuth=az, elevation=el)

    synth.SyntheticData(synth_field, radar)

if __name__ == '__main__':
    main()
