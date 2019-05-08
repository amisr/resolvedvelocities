# ResolveVectors.py

import tables
import numpy as np
import configparser

class ResolveVectors(object):
    def __init__(self):
        # read config file
        config = configparser.ConfigParser()
        config.read('config.ini')

        self.datafile = config.get('DEFAULT', 'DATAFILE')
        self.chirp = float(config.get('DEFAULT', 'CHIRP'))
        self.neMin = float(config.get('DEFAULT', 'NEMIN'))
        self.integration_time = float(config.get('DEFAULT', 'INTTIME'))

        # list of beam codes to use

    def read_data(self):
        # read data from standard AMISR fit files
        with tables.open_file(self.datafile,'r') as file:

            # time
            self.time = file.get_node('/Time/UnixTime')[:]

            # beam codes
            # define which beams to use (default is all)
            self.BeamCodes=file.get_node('/BeamCodes')[:,0]
            bm_idx = np.arange(0,len(self.BeamCodes))

            # geodetic location of each measurement
            self.alt = file.get_node('/Geomag/Altitude')[bm_idx,:].flatten()
            self.lat = file.get_node('/Geomag/Latitude')[bm_idx,:].flatten()
            self.lon = file.get_node('/Geomag/Longitude')[bm_idx,:].flatten()

            # geodetic k vectors
            self.ke = file.get_node('/Geomag/ke')[bm_idx,:].flatten()
            self.kn = file.get_node('/Geomag/kn')[bm_idx,:].flatten()
            self.kz = file.get_node('/Geomag/kz')[bm_idx,:].flatten()

            # line of sight velocity and error
            self.vlos = file.get_node('/FittedParams/Fits')[:,bm_idx,:,0,3].reshape((len(self.time[:,0]),len(self.alt)))
            # vlos1 = np.swapaxes(vlos1,0,1)
            self.dvlos = file.get_node('/FittedParams/Errors')[:,bm_idx,:,0,3].reshape((len(self.time[:,0]),len(self.alt)))
            # dvlos1 = np.swapaxes(dvlos1,0,1)

            # density (for filtering)
            self.ne = file.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))


    def filter_data(self):
        # filter and adjust data so it is appropriate for Bayesian reconstruction

        # add chirp to LoS velocity
        self.vlos = self.vlos + self.chirp

        # discard data with low density
        I = np.where((self.ne < self.neMin))
        self.vlos[I] = np.nan
        self.dvlos[I] = np.nan



    def transform(self):
        # transform k vectors from geodetic to geomagnetic
        # use apexpy
        pass

    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in some way by initial config file
        # each bin has a specified MLAT/MLON
        pass

    def compute_vectors(self):
        # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
        pass

    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        pass

    def save_output(self):
        # save output file
        pass