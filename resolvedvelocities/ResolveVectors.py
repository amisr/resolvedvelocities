# ResolveVectors.py

import tables
import numpy as np
import datetime as dt
import configparser
from apexpy import Apex


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

        # find where nans occur in input position arrays and remove them
        removed_nans = np.argwhere(np.isnan(self.alt)).flatten()
        glat = self.lat[np.isfinite(self.lat)]
        glon = self.lon[np.isfinite(self.lon)]
        galt = self.alt[np.isfinite(self.alt)]
        keg = self.ke[np.isfinite(self.ke)]
        kng = self.kn[np.isfinite(self.kn)]
        kzg = self.kz[np.isfinite(self.kz)]

        # intialize apex coordinates
        # time0 = dt.datetime.utcfromtimestamp(self.time[0,0])
        A = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))

        # find magnetic latitude and longitude
        mlat, mlon = A.geo2apex(glat, glon, galt/1000.)

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([keg, kng, kzg])
        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(glat, glon, galt/1000.)
        # find components of k for e1, e2, e3 basis vectors (Laundal and Richmond, 2016 eqn. 60)
        ke1 = np.einsum('ij,ij->j',kvec,d1)
        ke2 = np.einsum('ij,ij->j',kvec,d2)
        ke3 = np.einsum('ij,ij->j',kvec,d3)

        # reintroduce NANs
        # find indices where nans should be inserted in new arrays
        replace_nans = np.array([r-i for i,r in enumerate(removed_nans)])

        self.mlat = np.insert(mlat,replace_nans,np.nan)
        self.mlon = np.insert(mlon,replace_nans,np.nan)
        self.ke1 = np.insert(ke1,replace_nans,np.nan)
        self.ke2 = np.insert(ke2,replace_nans,np.nan)
        self.ke3 = np.insert(ke3,replace_nans,np.nan)



    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in some way by initial config file
        # each bin has a specified MLAT/MLON

        bin_edge_mlat = np.arange(65.0,69.0,0.5)
        self.data_bins = []
        for i in range(len(bin_edge_mlat)-1):
            center_mlat = (bin_edge_mlat[i]+bin_edge_mlat[i+1])/2.
            center_mlon = np.nanmean(self.mlon)
            idx = np.argwhere((self.mlat>=bin_edge_mlat[i]) & (self.mlat<bin_edge_mlat[i+1]))
            self.data_bins.append({'mlat':center_mlat,'mlon':center_mlon,'idx':idx.flatten()})


    # times: this is the the /Time/UnixTime array converted to datetime objects, so shape (num_recs,2)
    # integration_time:  this is the amount of seconds to integrate
    def get_integration_periods(self):
        self.integration_periods = []
        idx = []
        start_time = None
        num_times = len(self.time)
        for i,time_pair in enumerate(self.time):
            temp_start_time, temp_end_time = time_pair
            if start_time is None:
                start_time = temp_start_time
            time_diff = temp_end_time - start_time
            idx.append(i)

            if (time_diff >= self.integration_time) or (i == num_times -1):
                self.integration_periods.append({'start':start_time, 'end':temp_end_time, 'idx':np.array(idx)})
                idx = []
                start_time = None
                continue




    def compute_vectors(self):
        # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
        pass

    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        pass

    def save_output(self):
        # save output file
        pass