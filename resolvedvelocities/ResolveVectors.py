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
        self.chirp = eval(config.get('DEFAULT', 'CHIRP'))
        self.neMin = eval(config.get('DEFAULT', 'NEMIN'))
        self.integration_time = eval(config.get('DEFAULT', 'INTTIME'))
        self.covar = eval(config.get('DEFAULT', 'COVAR'))
        self.ppp = eval(config.get('DEFAULT', 'PPP'))
        self.minalt = eval(config.get('DEFAULT', 'MINALT'))
        self.maxalt = eval(config.get('DEFAULT', 'MAXALT'))
        self.minnumpoints = eval(config.get('DEFAULT', 'MINNUMPOINTS'))

        # self.covar=[1000.*1000.,1000.*1000.,5.*5.]
        # self.minalt = 150.
        # self.maxalt = 400.
        # self.ppp = [200.0,0.5,2000.0,100.0]
        # self.minnumpoints = 1

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

        # discard data outside of altitude range
        I = np.where(((self.alt < self.minalt*1000.) | (self.alt > self.maxalt*1000.)))
        self.vlos[I] = np.nan
        self.dvlos[I] = np.nan

        # discard data with "unexceptable" error
        #   - not sure about these conditions - they come from original vvels code but eliminate a lot of data points
        fracerrs = np.absolute(self.dvlos)/(np.absolute(self.vlos)+self.ppp[0])
        abserrs  = np.absolute(self.dvlos)
        I = np.where(((fracerrs > self.ppp[1]) & (abserrs > self.ppp[3])))
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
        self.Apex = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))

        # find magnetic latitude and longitude
        mlat, mlon = self.Apex.geo2apex(glat, glon, galt/1000.)

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([keg, kng, kzg])
        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, galt/1000.)
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
        
        Velocity = []
        VelocityCovariance = []

        # For each integration period and bin, calculate covarient components of drift velocity (Ve1, Ve2, Ve3)
        # loop over integration periods
        for ip in self.integration_periods:
            V_at_time = []
            SigV_at_time = []
            # loop over spatial bins
            for b in self.data_bins:

                # pull out the line of slight measurements for the time period and bins
                vlos = self.vlos[ip['idx'],b['idx'][:,np.newaxis]].flatten()
                dvlos = self.dvlos[ip['idx'],b['idx'][:,np.newaxis]].flatten()
                # pull out the k vectors for the bins and duplicate so they match the number of time measurements
                ke1 = np.repeat(self.ke1[b['idx']],len(ip['idx']))
                ke2 = np.repeat(self.ke2[b['idx']],len(ip['idx']))
                ke3 = np.repeat(self.ke3[b['idx']],len(ip['idx']))

                # remove nan data points
                finite = np.isfinite(vlos)
                vlos = vlos[finite]
                dvlos = dvlos[finite]
                ke1 = ke1[finite]
                ke2 = ke2[finite]
                ke3 = ke3[finite]


                A = np.array([ke1, ke2, ke3]).T
                SigmaE = np.diagflat(dvlos)
                SigmaV = np.diagflat(self.covar)

                try:
                    I = np.linalg.inv(np.einsum('jk,kl,ml->jm',A,SigmaV,A) + SigmaE)   # calculate I = (A*SigV*A.T + SigE)^-1
                    V = np.einsum('jk,lk,lm,m->j',SigmaV,A,I,vlos)      # calculate velocity estimate (Heinselman 2008 eqn 12)
                    SigV = np.linalg.inv(np.einsum('kj,kl,lm->jm',A,np.linalg.inv(SigmaE),A) + np.linalg.inv(SigmaV))       # calculate covariance of velocity estimate (Heinselman 2008 eqn 13)

                except np.linalg.LinAlgError:
                    V = np.full(3,np.nan)
                    SigV = np.full((3,3),np.nan)

                # if there are too few points for a valid reconstruction, set output to NAN
                if sum(finite) < self.minnumpoints:
                    V = np.full(3,np.nan)
                    SigV = np.full((3,3),np.nan)


                V_at_time.append(V)
                SigV_at_time.append(SigmaV)

            Velocity.append(V_at_time)
            VelocityCovariance.append(SigV_at_time)

        self.Velocity = np.array(Velocity)
        self.VelocityCovariance = np.array(VelocityCovariance)


        # calculate electric field
        bin_mlat = [b['mlat'] for b in self.data_bins]
        bin_mlon = [b['mlon'] for b in self.data_bins]

        # find Be3 value at each output bin location
        Be3, __, __, __ = self.Apex.bvectors_apex(bin_mlat,bin_mlon,300.,coords='apex')
        # Be3 = np.full(plat_out1.shape,1.0)        # set Be3 array to 1.0 - useful for debugging linear algebra

        # form rotation array
        R = np.einsum('i,jk->ijk',Be3,np.array([[0,-1,0],[1,0,0],[0,0,0]]))
        # Calculate contravarient components of electric field (Ed1, Ed2, Ed3)
        self.ElectricField = np.einsum('ijk,...ik->...ij',R,self.Velocity)
        # Calculate electric field covariance matrix (SigE = R*SigV*R.T)
        self.ElectricFieldCovariance = np.einsum('ijk,...ikl,iml->...ijm',R,self.VelocityCovariance,R)

        # print self.Velocity.shape, self.ElectricField.shape


    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        pass

    def save_output(self):
        # save output file
        pass