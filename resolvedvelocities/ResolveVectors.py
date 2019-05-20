# ResolveVectors.py

import tables
import numpy as np
import datetime as dt
import configparser
from apexpy import Apex

# TODO:
# - Make output bin coordinates independent arrays
# - Use consistent notation for input/output parameters

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

        # remove nans from input position arrays (apexpy can't handle input nans)
        glat = self.lat[np.isfinite(self.lat)]
        glon = self.lon[np.isfinite(self.lon)]
        galt = self.alt[np.isfinite(self.alt)]
        keg = self.ke[np.isfinite(self.ke)]
        kng = self.kn[np.isfinite(self.kn)]
        kzg = self.kz[np.isfinite(self.kz)]

        # intialize apex coordinates
        self.Apex = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))

        # find magnetic latitude and longitude
        mlat, mlon = self.Apex.geo2apex(glat, glon, galt/1000.)

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([keg, kng, kzg]).T
        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, galt/1000.)
        d = np.array([d1,d2,d3]).T
        # find components of k for e1, e2, e3 basis vectors (Laundal and Richmond, 2016 eqn. 60)
        A = np.einsum('ij,ijk->ik', kvec, d)

        self.mlat = mlat
        self.mlon = mlon
        self.A = A




    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in some way by initial config file
        # each bin has a specified MLAT/MLON

        bin_edge_mlat = np.arange(64.0,68.0,0.25)
        self.bin_mlat = (bin_edge_mlat[:-1]+bin_edge_mlat[1:])/2.
        self.bin_mlon = np.full(self.bin_mlat.shape, np.nanmean(self.mlon))
        self.bin_idx = [np.argwhere((self.mlat>=bin_edge_mlat[i]) & (self.mlat<bin_edge_mlat[i+1])).flatten() for i in range(len(bin_edge_mlat)-1)]



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
            Vel = []
            SigmaV = []
            # loop over spatial bins
            for bidx in self.bin_idx:

                # pull out the line of slight measurements for the time period and bins
                vlos = self.vlos[ip['idx'],bidx[:,np.newaxis]].flatten()
                dvlos = self.dvlos[ip['idx'],bidx[:,np.newaxis]].flatten()
                # pull out the k vectors for the bins and duplicate so they match the number of time measurements
                A = np.repeat(self.A[bidx], len(ip['idx']), axis=0)

                # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
                V, SigV = vvels(vlos, dvlos, A, self.covar, minnumpoints=self.minnumpoints)

                # append vector and coviarience matrix
                Vel.append(V)
                SigmaV.append(SigV)

            Velocity.append(Vel)
            VelocityCovariance.append(SigmaV)

        self.Velocity = np.array(Velocity)
        self.VelocityCovariance = np.array(VelocityCovariance)


        # calculate electric field
        # find Be3 value at each output bin location
        Be3, __, __, __ = self.Apex.bvectors_apex(self.bin_mlat,self.bin_mlon,200.,coords='apex')
        # Be3 = np.full(plat_out1.shape,1.0)        # set Be3 array to 1.0 - useful for debugging linear algebra

        # form rotation array
        R = np.einsum('i,jk->ijk',Be3,np.array([[0,-1,0],[1,0,0],[0,0,0]]))
        # Calculate contravarient components of electric field (Ed1, Ed2, Ed3)
        self.ElectricField = np.einsum('ijk,...ik->...ij',R,self.Velocity)
        # Calculate electric field covariance matrix (SigE = R*SigV*R.T)
        self.ElectricFieldCovariance = np.einsum('ijk,...ikl,iml->...ijm',R,self.VelocityCovariance,R)



    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes

        alt = 200.  # for now, just calculate vectors at a set altitude

        # calculate bin locations in geodetic coordinates
        gdlat, gdlon, err = self.Apex.apex2geo(self.bin_mlat, self.bin_mlon, alt)
        self.gdlat = gdlat
        self.gdlon = gdlon
        self.gdalt = np.full(gdlat.shape, alt)

        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(self.bin_mlat, self.bin_mlon, alt, coords='apex')

        e = np.array([e1,e2,e3]).T
        self.Velocity_gd = np.einsum('ijk,...ik->...ij',e,self.Velocity)
        self.VelocityCovariance_gd = np.einsum('ijk,...ikl,iml->...ijm',e,self.VelocityCovariance,e)

        d = np.array([d1,d2,d3]).T
        self.ElectricField_gd = np.einsum('ijk,...ik->...ij',d,self.ElectricField)
        self.ElectricFieldCovariance_gd = np.einsum('ijk,...ikl,iml->...ijm',d,self.ElectricFieldCovariance,d)

        # calculate vector magnitude
        self.Vgd_mag = np.linalg.norm(self.Velocity_gd,axis=-1)
        self.Egd_mag = np.linalg.norm(self.ElectricField_gd,axis=-1)


    def save_output(self):

        out_time = [[t['start'],t['end']] for t in self.integration_periods]

        # save output file
        filename = 'test_vvels.h5'

        with tables.open_file(filename,mode='w') as file:
            file.create_array('/','UnixTime',out_time)
            file.set_node_attr('/UnixTime', 'TITLE', 'Unix Time')
            file.set_node_attr('/UnixTime', 'Size', 'Nrecords x 2 (start and end of integration)')
            file.set_node_attr('/UnixTime', 'Units', 's')
            file.create_group('/', 'Magnetic')
            file.create_group('/', 'Geographic')
            file.create_array('/Magnetic', 'MagneticLatitude', self.bin_mlat)
            file.set_node_attr('/Magnetic/MagneticLatitude', 'TITLE', 'Magnetic Latitude')
            file.set_node_attr('/Magnetic/MagneticLatitude', 'Size', 'Nbins')
            file.create_array('/Magnetic','MagneticLongitude', self.bin_mlon)
            file.set_node_attr('/Magnetic/MagneticLongitude', 'TITLE', 'Magnetic Longitude')
            file.set_node_attr('/Magnetic/MagneticLongitude', 'Size', 'Nbins')
            file.create_array('/Magnetic', 'Velocity', self.Velocity)
            file.set_node_attr('/Magnetic/Velocity', 'TITLE', 'Plama Drift Velocity')
            file.set_node_attr('/Magnetic/Velocity', 'Size', 'Nrecords x Nbins x 3 (Ve1, Ve2, Ve3)')
            file.set_node_attr('/Magnetic/Velocity', 'Units', 'm/s')
            file.create_array('/Magnetic','SigmaV', self.VelocityCovariance)
            file.set_node_attr('/Magnetic/SigmaV', 'TITLE', 'Velocity Covariance Matrix')
            file.set_node_attr('/Magnetic/SigmaV', 'Size', 'Nrecords x Nbins x 3 x 3')
            file.set_node_attr('/Magnetic/SigmaV', 'Units', 'm/s')
            file.create_array('/Magnetic','ElectricField',self.ElectricField)
            file.set_node_attr('/Magnetic/ElectricField', 'TITLE', 'Convection Electric Field')
            file.set_node_attr('/Magnetic/ElectricField', 'Size', 'Nrecords x Nbins x 3 (Ed1, Ed2, Ed3)')
            file.set_node_attr('/Magnetic/ElectricField', 'Units', 'V/m')
            file.create_array('/Magnetic','SigmaE',self.ElectricFieldCovariance)
            file.set_node_attr('/Magnetic/SigmaE', 'TITLE', 'Electric Field Covariance Matrix')
            file.set_node_attr('/Magnetic/SigmaE', 'Size', 'Nrecords x Nbins x 3 x 3')
            file.set_node_attr('/Magnetic/SigmaE', 'Units', 'V/m')


def vvels(vlos, dvlos, A, cov, minnumpoints=1):
    # implimentation of Heinselman and Nicolls 2008 vector velocity algorithm

    # remove nan data points
    finite = np.isfinite(vlos)
    vlos = vlos[finite]
    dvlos = dvlos[finite]
    A = A[finite]

    SigmaE = np.diagflat(dvlos**2)
    SigmaV = np.diagflat(cov)

    try:
        I = np.linalg.inv(np.einsum('jk,kl,ml->jm',A,SigmaV,A) + SigmaE)   # calculate I = (A*SigV*A.T + SigE)^-1
        V = np.einsum('jk,lk,lm,m->j',SigmaV,A,I,vlos)      # calculate velocity estimate (Heinselman 2008 eqn 12)
        SigV = np.linalg.inv(np.einsum('kj,kl,lm->jm',A,np.linalg.inv(SigmaE),A) + np.linalg.inv(SigmaV))       # calculate covariance of velocity estimate (Heinselman 2008 eqn 13)

    except np.linalg.LinAlgError:
        V = np.full(3,np.nan)
        SigV = np.full((3,3),np.nan)

    # if there are too few points for a valid reconstruction, set output to NAN
    if sum(finite) < minnumpoints:
        V = np.full(3,np.nan)
        SigV = np.full((3,3),np.nan)

    return V, SigV