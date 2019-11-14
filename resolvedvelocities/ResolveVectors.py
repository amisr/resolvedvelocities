# ResolveVectors.py

import tables
import numpy as np
import datetime as dt
import os
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import platform
import getpass
import socket
from apexpy import Apex
from .marp import Marp
from scipy.spatial import Delaunay


class ResolveVectors(object):
    def __init__(self, config):
        self.configfile = config
        self.read_config(self.configfile)

    def read_config(self, config_file):
        # read config file
        config = configparser.ConfigParser()
        config.read(config_file)

        self.datafile = config.get('DEFAULT', 'DATAFILE')
        self.outfilename = config.get('DEFAULT', 'OUTFILENAME')
        self.chirp = eval(config.get('DEFAULT', 'CHIRP'))
        self.covar = eval(config.get('DEFAULT', 'COVAR'))
        self.altlim = eval(config.get('DEFAULT', 'ALTLIM'))
        self.nelim = eval(config.get('DEFAULT', 'NELIM'))
        self.chi2lim = eval(config.get('DEFAULT', 'CHI2LIM'))
        self.goodfitcode = eval(config.get('DEFAULT', 'GOODFITCODE'))
        self.binvert = eval(config.get('DEFAULT', 'BINVERT'))
        self.outalt = eval(config.get('DEFAULT', 'OUTALT'))
        self.minnumpoints = eval(config.get('DEFAULT', 'MINNUMPOINTS'))

        if config.has_option('DEFAULT', 'UPB_BEAMCODE'):
            self.upB_beamcode = config.getint('DEFAULT', 'UPB_BEAMCODE')
        else:
            self.upB_beamcode = None
        if config.has_option('DEFAULT', 'IONUP'):
            self.ionup = config.get('DEFAULT', 'IONUP')
        else:
            self.ionup = None
        if config.has_option('DEFAULT', 'USE_BEAMS'):
            self.use_beams = eval(config.get('DEFAULT', 'USE_BEAMS'))
        else:
            self.use_beams = None
        if config.has_option('DEFAULT', 'INTTIME'):
            self.integration_time = config.getfloat('DEFAULT', 'INTTIME')
        else:
            self.integration_time = None


    def read_data(self):
        # read data from standard AMISR fit files
        with tables.open_file(self.datafile,'r') as infile:

            # time
            self.time = infile.get_node('/Time/UnixTime')[:]

            # site
            lat = infile.get_node('/Site/Latitude').read()
            lon = infile.get_node('/Site/Longitude').read()
            alt = infile.get_node('/Site/Altitude').read()
            self.site = np.array([lat, lon, alt/1000.])

            # define which beams to use (default is all)
            self.BeamCodes=infile.get_node('/BeamCodes')[:,0]
            if self.use_beams:
                bm_idx = np.array([i for i,b in enumerate(self.BeamCodes) if b in self.use_beams])
            else:
                bm_idx = np.arange(0,len(self.BeamCodes))

            # geodetic location of each measurement
            self.alt = infile.get_node('/Geomag/Altitude')[bm_idx,:].flatten()
            self.lat = infile.get_node('/Geomag/Latitude')[bm_idx,:].flatten()
            self.lon = infile.get_node('/Geomag/Longitude')[bm_idx,:].flatten()

            # geodetic k vectors
            self.ke = infile.get_node('/Geomag/ke')[bm_idx,:].flatten()
            self.kn = infile.get_node('/Geomag/kn')[bm_idx,:].flatten()
            self.kz = infile.get_node('/Geomag/kz')[bm_idx,:].flatten()

            # line of sight velocity and error
            self.vlos = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,0,3].reshape((len(self.time[:,0]),len(self.alt)))
            self.dvlos = infile.get_node('/FittedParams/Errors')[:,bm_idx,:,0,3].reshape((len(self.time[:,0]),len(self.alt)))

            # chi2 and fitcode (for filtering poor quality data)
            self.chi2 = infile.get_node('/FittedParams/FitInfo/chi2')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))
            self.fitcode = infile.get_node('/FittedParams/FitInfo/fitcode')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))

            # density (for filtering and ion upflow)
            self.ne = infile.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))

            # temperature (for ion upflow)
            self.Te = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,5,1].reshape((len(self.time[:,0]),len(self.alt)))
            Ts = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:5,1]
            frac = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:5,0]
            self.Ti = np.sum(Ts*frac,axis=-1).reshape((len(self.time[:,0]),len(self.alt)))

            # get up-B beam velocities for ion outflow correction
            if self.upB_beamcode:
                upB_idx = np.argwhere(self.BeamCodes==self.upB_beamcode).flatten()
                upB_alt = infile.get_node('/Geomag/Altitude')[upB_idx,:].flatten()
                upB_vlos = infile.get_node('/FittedParams/Fits')[:,upB_idx,:,0,3].reshape((len(self.time[:,0]),len(upB_alt)))
                upB_dvlos = infile.get_node('/FittedParams/Errors')[:,upB_idx,:,0,3].reshape((len(self.time[:,0]),len(upB_alt)))
                self.upB = {'alt':upB_alt, 'vlos':upB_vlos, 'dvlos':upB_dvlos}



    def filter_data(self):
        # filter and adjust data so it is appropriate for Bayesian reconstruction

        # add chirp to LoS velocity
        self.vlos = self.vlos + self.chirp

        # discard data with low density
        I = np.where((self.ne < self.nelim[0]) | (self.ne > self.nelim[1]))
        self.vlos[I] = np.nan
        self.dvlos[I] = np.nan

        # discard data outside of altitude range
        I = np.where((self.alt < self.altlim[0]*1000.) | (self.alt > self.altlim[1]*1000.))
        self.vlos[:,I] = np.nan
        self.dvlos[:,I] = np.nan

        # discard data with extremely high or extremely low chi2 values
        I = np.where((self.chi2 < self.chi2lim[0]) | (self.chi2 > self.chi2lim[1]))
        self.vlos[I] = np.nan
        self.dvlos[I] = np.nan

        # discard data with poor fitcode (fitcodes 1-4 denote solution found, anything else should not be used)
        I = np.where(~np.isin(self.fitcode, self.goodfitcode))
        self.vlos[I] = np.nan
        self.vlos[I] = np.nan


    def transform(self):
        # transform k vectors from geodetic to geomagnetic

        # find indices where nans will be removed and should be inserted in new arrays
        replace_nans = np.array([r-i for i,r in enumerate(np.argwhere(np.isnan(self.alt)).flatten())])

        glat = self.lat[np.isfinite(self.lat)]
        glon = self.lon[np.isfinite(self.lon)]
        galt = self.alt[np.isfinite(self.alt)]/1000.

        A = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))
        alat, alon = A.geo2apex(glat, glon, galt)
        print(np.mean(alat), np.mean(alon))

        # intialize apex coordinates
        # self.Apex = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))
        # self.marp = Marp(date=dt.datetime.utcfromtimestamp(self.time[0,0]), lam0=0., phi0=0.)
        self.marp = Marp(date=dt.datetime.utcfromtimestamp(self.time[0,0]), lam0=np.mean(alat), phi0=np.mean(alon))

        # find magnetic latitude and longitude
        # mlat, mlon = self.Apex.geo2apex(glat, glon, galt)
        mlat, mlon = self.marp.geo2marp(glat, glon, galt)
        self.mlat = np.insert(mlat,replace_nans,np.nan)
        self.mlon = np.insert(mlon,replace_nans,np.nan)
        # print(self.mlat, self.mlon)

        # apex basis vectors in geodetic coordinates [e n u]
        # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, galt)
        d1, d2, d3, e1, e2, e3 = self.marp.basevectors_marp(glat, glon, galt)
        d1 = np.insert(d1,replace_nans,np.nan,axis=1)
        d2 = np.insert(d2,replace_nans,np.nan,axis=1)
        d3 = np.insert(d3,replace_nans,np.nan,axis=1)
        e1 = np.insert(e1,replace_nans,np.nan,axis=1)
        e2 = np.insert(e2,replace_nans,np.nan,axis=1)
        e3 = np.insert(e3,replace_nans,np.nan,axis=1)
        e = np.array([e1,e2,e3]).T

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([self.ke, self.kn, self.kz]).T

        # find components of k for d1, d2, d3 base vectors (Laundal and Richmond, 2016 eqn. 60)
        self.A = np.einsum('ij,ijk->ik', kvec, e)

        # calculate scaling factor D, used for ion outflow correction (Richmond, 1995 eqn. 3.15)
        d1_cross_d2 = np.cross(d1.T,d2.T).T
        self.D = np.sqrt(np.sum(d1_cross_d2**2,axis=0))



    def ion_upflow_correction(self):

        # correct the los velocities for the entire array at each time
        for t in range(len(self.time)):
            if not self.ionup:
                continue
            elif self.ionup == 'UPB':
                # interpolate velocities from up B beam to all other measurements
                vion, dvion = lin_interp(self.alt, self.upB['alt'], self.upB['vlos'][t], self.upB['dvlos'][t])
            elif self.ionup == 'EMP':
                # use empirical method to find ion upflow
                # NOTE: NOT DEVELOPED YET!!!
                vion, dvion = ion_upflow(self.Te, self.Ti, self.ne, self.alt)

            # LoS velocity correction to remove ion upflow
            self.vlos[t] = self.vlos[t] + self.A[:,2]/self.D*vion
            # corrected error in new LoS velocities
            self.dvlos[t] = np.sqrt(self.dvlos[t]**2 + self.A[:,2]**2/self.D**2*dvion**2)



    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in config file by a list of bin vericies in apex magnetic coordinates
        # the center of ecah bin is defined as the average of the verticies

        self.bin_mlat = []
        self.bin_mlon = []
        self.bin_idx = []
        for vert in self.binvert:
            vert = np.array(vert)
            hull = Delaunay(vert)

            self.bin_mlat.append(np.nanmean(vert[:,0]))
            self.bin_mlon.append(np.nanmean(vert[:,1]))
            self.bin_idx.append(np.argwhere(hull.find_simplex(np.array([self.mlat, self.mlon]).T)>=0).flatten())



    def get_integration_periods(self):

        if not self.integration_time:
            # if no integration time specified, use original time periods of input files
            self.int_period = self.time
            self.int_idx = range(len(self.time))

        else:
            # if an integration time is given, calculate new time periods
            self.int_period = []
            self.int_idx = []

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
                    self.int_period.append([start_time, temp_end_time])
                    self.int_idx.append(np.array(idx))
                    idx = []
                    start_time = None
                    continue

            self.int_period = np.array(self.int_period)


    def compute_vector_velocity(self):
        # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors

        Velocity = []
        VelocityCovariance = []
        ChiSquared = []

        # For each integration period and bin, calculate covarient components of drift velocity (Ve1, Ve2, Ve3)
        # loop over integration periods
        for tidx in self.int_idx[:1]:
            Vel = []
            SigmaV = []
            Chi2 = []
            # loop over spatial bins
            for bidx in self.bin_idx:

                # pull out the line of slight measurements for the time period and bins
                vlos = self.vlos[tidx,bidx[:,np.newaxis]].flatten()
                dvlos = self.dvlos[tidx,bidx[:,np.newaxis]].flatten()
                # pull out the k vectors for the bins and duplicate so they match the number of time measurements
                if self.integration_time:
                    A = np.repeat(self.A[bidx], len(tidx), axis=0)
                else:
                    # if no post integraiton, k vectors do not need to be duplicated
                    A = self.A[bidx]

                # print(A)

                # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
                V, SigV, chi2 = vvels(vlos, dvlos, A, self.covar, minnumpoints=self.minnumpoints)

                # append vector and coviarience matrix
                Vel.append(V)
                SigmaV.append(SigV)
                Chi2.append(chi2)

            Velocity.append(Vel)
            VelocityCovariance.append(SigmaV)
            ChiSquared.append(Chi2)

        self.Velocity = np.array(Velocity)
        self.VelocityCovariance = np.array(VelocityCovariance)
        self.ChiSquared = np.array(ChiSquared)


    def compute_electric_field(self):
        # calculate electric field

        # find Be3 value at each output bin location
        # NOTE: Be3 is constant along magnetic field lines, so the altitude chosen here doesn't matter
        # Be3, __, __, __ = self.Apex.bvectors_apex(self.bin_mlat,self.bin_mlon,300.,coords='apex')
        Be3 = np.full(len(self.bin_mlat),1.0)        # set Be3 array to 1.0 - useful for debugging linear algebra

        # form rotation array
        R = np.einsum('i,jk->ijk',Be3,np.array([[0,-1,0],[1,0,0],[0,0,0]]))
        # Calculate contravarient components of electric field (Ed1, Ed2, Ed3)
        self.ElectricField = np.einsum('ijk,...ik->...ij',R,self.Velocity)
        # Calculate electric field covariance matrix (SigE = R*SigV*R.T)
        self.ElectricFieldCovariance = np.einsum('ijk,...ikl,iml->...ijm',R,self.VelocityCovariance,R)



    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        # altitudes are defined by config file

        hbins = len(self.bin_mlat)
        vbins = len(self.outalt)
        mlat = np.tile(self.bin_mlat,vbins)
        mlon = np.tile(self.bin_mlon,vbins)
        alt = np.repeat(self.outalt, hbins)

        # calculate bin locations in geodetic coordinates
        # glat, glon, err = self.Apex.apex2geo(mlat, mlon, alt)
        glat, glon, err = self.marp.marp2geo(mlat, mlon, alt)
        self.bin_glat = glat.reshape((vbins,hbins))
        self.bin_glon = glon.reshape((vbins,hbins))
        self.bin_galt = alt.reshape((vbins,hbins))

        # apex basis vectors in geodetic coordinates [e n u]
        # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(mlat, mlon, alt, coords='apex')
        # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, alt)
        d1, d2, d3, e1, e2, e3 = self.marp.basevectors_marp(glat, glon, alt)

        # Ve3 and Ed3 should be 0 because VE and E should not have components parallel to B.
        # To force this, set e3 = 0 and d3 = 0
        # Because E was derived by taking the cross product between B and V, setting d3=0 SHOULD be redundant
        e3 = np.zeros(e3.shape)
        d3 = np.zeros(d3.shape)

        # caluclate plasma drift velocity geodetic components
        e = np.array([e1,e2,e3]).T.reshape((vbins,hbins,3,3))
        self.Velocity_gd = np.einsum('hijk,...ik->...hij',e,self.Velocity)
        self.VelocityCovariance_gd = np.einsum('hijk,...ikl,himl->...hijm',e,self.VelocityCovariance,e)

        # calculate electric field geodetic components
        d = np.array([d1,d2,d3]).T.reshape((vbins,hbins,3,3))
        self.ElectricField_gd = np.einsum('hijk,...ik->...hij',d,self.ElectricField)
        self.ElectricFieldCovariance_gd = np.einsum('hijk,...ikl,himl->...hijm',d,self.ElectricFieldCovariance,d)

        # calculate vector magnitude and direction
        north = -e2.T.reshape((vbins,hbins,3))
        self.Vgd_mag, self.Vgd_mag_err, self.Vgd_dir, self.Vgd_dir_err = magnitude_direction(self.Velocity_gd, self.VelocityCovariance_gd, north)
        self.Egd_mag, self.Egd_mag_err, self.Egd_dir, self.Egd_dir_err = magnitude_direction(self.ElectricField_gd, self.ElectricFieldCovariance_gd, north)


    def save_output(self):
        # TODO: come up with a better way to manage all this

        # save output file
        with tables.open_file(self.outfilename,mode='w') as outfile:

            # copy some groups directly from fitted input file
            with tables.open_file(self.datafile, mode='r') as infile:
                outfile.copy_children(infile.get_node('/Site'), outfile.create_group('/','Site'))
                if not self.integration_time:
                    outfile.copy_children(infile.get_node('/Time'), outfile.create_group('/','Time'))
                else:
                    outfile.create_group('/','Time')
                    year, month, day, doy, dtime, mlt = self.create_time_arrays()

                    outfile.create_array('/Time', 'UnixTime', self.int_period)
                    outfile.set_node_attr('/Time/UnixTime', 'TITLE', 'UnixTime')
                    outfile.set_node_attr('/Time/UnixTime', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/UnixTime', 'Unit', 'Seconds')

                    outfile.create_array('/Time', 'Year', year)
                    outfile.set_node_attr('/Time/Year', 'TITLE', 'Year')
                    outfile.set_node_attr('/Time/Year', 'Size', 'Nrecords x 2 (Start and end of integration')

                    outfile.create_array('/Time', 'Month', month)
                    outfile.set_node_attr('/Time/Month', 'TITLE', 'Month')
                    outfile.set_node_attr('/Time/Month', 'Size', 'Nrecords x 2 (Start and end of integration')

                    outfile.create_array('/Time', 'Day', day)
                    outfile.set_node_attr('/Time/Day', 'TITLE', 'Day of Month')
                    outfile.set_node_attr('/Time/Day', 'Size', 'Nrecords x 2 (Start and end of integration')

                    outfile.create_array('/Time', 'doy', doy)
                    outfile.set_node_attr('/Time/doy', 'TITLE', 'Day of Year')
                    outfile.set_node_attr('/Time/doy', 'Size', 'Nrecords x 2 (Start and end of integration')

                    outfile.create_array('/Time', 'dtime', dtime)
                    outfile.set_node_attr('/Time/dtime', 'TITLE', 'Decimal Time')
                    outfile.set_node_attr('/Time/dtime', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/dtime', 'Unit', 'UT Hour')

                    outfile.create_array('/Time', 'MagneticLocalTimeSite', mlt)
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'TITLE', 'Magnetic Local Time')
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'Unit', 'UT Hour')

            outfile.create_group('/', 'Magnetic')

            outfile.create_array('/Magnetic', 'MagneticLatitude', self.bin_mlat)
            outfile.set_node_attr('/Magnetic/MagneticLatitude', 'TITLE', 'Magnetic Latitude')
            outfile.set_node_attr('/Magnetic/MagneticLatitude', 'Size', 'Nbins')

            outfile.create_array('/Magnetic','MagneticLongitude', self.bin_mlon)
            outfile.set_node_attr('/Magnetic/MagneticLongitude', 'TITLE', 'Magnetic Longitude')
            outfile.set_node_attr('/Magnetic/MagneticLongitude', 'Size', 'Nbins')

            outfile.create_array('/Magnetic', 'Velocity', self.Velocity)
            outfile.set_node_attr('/Magnetic/Velocity', 'TITLE', 'Plama Drift Velocity')
            outfile.set_node_attr('/Magnetic/Velocity', 'Size', 'Nrecords x Nbins x 3 (Ve1, Ve2, Ve3)')
            outfile.set_node_attr('/Magnetic/Velocity', 'Units', 'm/s')

            outfile.create_array('/Magnetic','SigmaV', self.VelocityCovariance)
            outfile.set_node_attr('/Magnetic/SigmaV', 'TITLE', 'Velocity Covariance Matrix')
            outfile.set_node_attr('/Magnetic/SigmaV', 'Size', 'Nrecords x Nbins x 3 x 3')
            outfile.set_node_attr('/Magnetic/SigmaV', 'Units', 'm/s')

            outfile.create_array('/Magnetic','ElectricField',self.ElectricField)
            outfile.set_node_attr('/Magnetic/ElectricField', 'TITLE', 'Convection Electric Field')
            outfile.set_node_attr('/Magnetic/ElectricField', 'Size', 'Nrecords x Nbins x 3 (Ed1, Ed2, Ed3)')
            outfile.set_node_attr('/Magnetic/ElectricField', 'Units', 'V/m')

            outfile.create_array('/Magnetic','SigmaE',self.ElectricFieldCovariance)
            outfile.set_node_attr('/Magnetic/SigmaE', 'TITLE', 'Electric Field Covariance Matrix')
            outfile.set_node_attr('/Magnetic/SigmaE', 'Size', 'Nrecords x Nbins x 3 x 3')
            outfile.set_node_attr('/Magnetic/SigmaE', 'Units', 'V/m')

            outfile.create_array('/Magnetic', 'Chi2', self.ChiSquared)
            outfile.set_node_attr('/Magnetic/Chi2', 'TITLE', 'Reduced Chi-Squared')
            outfile.set_node_attr('/Magnetic/Chi2', 'Size', 'Nrecords x Nbins')

            outfile.create_group('/', 'Geographic')

            outfile.create_array('/Geographic', 'GeographicLatitude', self.bin_glat)
            outfile.set_node_attr('/Geographic/GeographicLatitude', 'TITLE', 'Geographic Latitude')
            outfile.set_node_attr('/Geographic/GeographicLatitude', 'Size', 'Nalt x Nbins')

            outfile.create_array('/Geographic','GeographicLongitude', self.bin_glon)
            outfile.set_node_attr('/Geographic/GeographicLongitude', 'TITLE', 'Geographic Longitude')
            outfile.set_node_attr('/Geographic/GeographicLongitude', 'Size', 'Nalt x Nbins')

            outfile.create_array('/Geographic','GeographicAltitude', self.bin_galt)
            outfile.set_node_attr('/Geographic/GeographicAltitude', 'TITLE', 'Geographic Altitude')
            outfile.set_node_attr('/Geographic/GeographicAltitude', 'Size', 'Nalt x Nbins')
            outfile.set_node_attr('/Geographic/GeographicAltitude', 'Units', 'km')

            outfile.create_array('/Geographic', 'Velocity', self.Velocity_gd)
            outfile.set_node_attr('/Geographic/Velocity', 'TITLE', 'Plama Drift Velocity')
            outfile.set_node_attr('/Geographic/Velocity', 'Size', 'Nrecords x Nalt x Nbins x 3 (East, North, Up)')
            outfile.set_node_attr('/Geographic/Velocity', 'Units', 'm/s')

            outfile.create_array('/Geographic','SigmaV', self.VelocityCovariance_gd)
            outfile.set_node_attr('/Geographic/SigmaV', 'TITLE', 'Velocity Covariance Matrix')
            outfile.set_node_attr('/Geographic/SigmaV', 'Size', 'Nrecords x Nalt x Nbins x 3 x 3')
            outfile.set_node_attr('/Geographic/SigmaV', 'Units', 'm/s')

            outfile.create_array('/Geographic','Vmag',self.Vgd_mag)
            outfile.set_node_attr('/Geographic/Vmag', 'TITLE', 'Velocity Magnitude')
            outfile.set_node_attr('/Geographic/Vmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/Vmag', 'Units', 'm/s')

            outfile.create_array('/Geographic','errVmag',self.Vgd_mag_err)
            outfile.set_node_attr('/Geographic/errVmag', 'TITLE', 'Velocity Magnitude Error')
            outfile.set_node_attr('/Geographic/errVmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/errVmag', 'Units', 'm/s')

            outfile.create_array('/Geographic','Vdir',self.Vgd_dir)
            outfile.set_node_attr('/Geographic/Vdir', 'TITLE', 'Velocity Direction Angle East of North Magnetic Meridian (-e2)')
            outfile.set_node_attr('/Geographic/Vdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/Vdir', 'Units', 'Degrees')

            outfile.create_array('/Geographic','errVdir',self.Vgd_dir_err)
            outfile.set_node_attr('/Geographic/errVdir', 'TITLE', 'Error in Velocity Direction')
            outfile.set_node_attr('/Geographic/errVdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/errVdir', 'Units', 'Degrees')

            outfile.create_array('/Geographic','ElectricField',self.ElectricField_gd)
            outfile.set_node_attr('/Geographic/ElectricField', 'TITLE', 'Convection Electric Field')
            outfile.set_node_attr('/Geographic/ElectricField', 'Size', 'Nrecords x Nalt x Nbins x 3 (East, North, Up)')
            outfile.set_node_attr('/Geographic/ElectricField', 'Units', 'V/m')

            outfile.create_array('/Geographic','SigmaE',self.ElectricFieldCovariance_gd)
            outfile.set_node_attr('/Geographic/SigmaE', 'TITLE', 'Electric Field Covariance Matrix')
            outfile.set_node_attr('/Geographic/SigmaE', 'Size', 'Nrecords x Nalt x Nbins x 3 x 3')
            outfile.set_node_attr('/Geographic/SigmaE', 'Units', 'V/m')

            outfile.create_array('/Geographic','Emag',self.Egd_mag)
            outfile.set_node_attr('/Geographic/Emag', 'TITLE', 'Electric Field Magnitude')
            outfile.set_node_attr('/Geographic/Emag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/Emag', 'Units', 'V/m')

            outfile.create_array('/Geographic','errEmag',self.Egd_mag_err)
            outfile.set_node_attr('/Geographic/errEmag', 'TITLE', 'Electric Field Magnitude Error')
            outfile.set_node_attr('/Geographic/errEmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/errEmag', 'Units', 'V/m')

            outfile.create_array('/Geographic','Edir',self.Egd_dir)
            outfile.set_node_attr('/Geographic/Edir', 'TITLE', 'Electric Field Direction Angle East of North Magnetic Meridian (-e2)')
            outfile.set_node_attr('/Geographic/Edir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/Edir', 'Units', 'Degrees')

            outfile.create_array('/Geographic','errEdir',self.Egd_dir_err)
            outfile.set_node_attr('/Geographic/errEdir', 'TITLE', 'Error in Electric Field Direction')
            outfile.set_node_attr('/Geographic/errEdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geographic/errEdir', 'Units', 'Degrees')

            outfile.create_group('/', 'ProcessingParams')

            # outfile.create_array('/ProcessingParams', 'ApexYear', self.Apex.year)
            outfile.create_array('/ProcessingParams', 'ApexYear', self.marp.year)
            outfile.set_node_attr('/ProcessingParams/ApexYear', 'TITLE', 'Decimal Year used for IGRF Model')

            # outfile.create_array('/ProcessingParams', 'ApexRefHeight', self.Apex.refh)
            outfile.create_array('/ProcessingParams', 'ApexRefHeight', self.marp.refh)
            outfile.set_node_attr('/ProcessingParams/ApexRefHeight', 'TITLE', 'Reference height used for Apex coordinates')
            outfile.set_node_attr('/ProcessingParams/ApexRefHeight', 'Units', 'km')


            outfile.create_array('/ProcessingParams', 'InputFile', self.datafile.encode('utf-8'))

            # Save computer information and config file
            # Computer information:
            PythonVersion   = platform.python_version()
            Type            = platform.machine()
            System          = '{} {} {}'.format(platform.system(),platform.release(),platform.version())
            User            = getpass.getuser()
            Hostname        = platform.node()
            if len(Hostname) == 0:
                Hostname = socket.gethostname()

            outfile.create_group('/ProcessingParams', 'ComputerInfo')

            outfile.create_array('/ProcessingParams/ComputerInfo', 'PythonVersion', PythonVersion.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ComputerInfo', 'Type', Type.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ComputerInfo', 'System', System.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ComputerInfo', 'User', User.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ComputerInfo', 'Host', Hostname.encode('utf-8'))

            Path = os.path.dirname(os.path.abspath(self.configfile))
            Name = os.path.basename(self.configfile)
            with open(self.configfile, 'r') as f:
                Contents = ''.join(f.readlines())

            outfile.create_group('/ProcessingParams', 'ConfigFile')

            outfile.create_array('/ProcessingParams/ConfigFile', 'Name', Name.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ConfigFile', 'Path', Path.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ConfigFile', 'Contents', Contents.encode('utf-8'))


    def create_time_arrays(self):
        time_array = np.array([[dt.datetime.utcfromtimestamp(t[0]), dt.datetime.utcfromtimestamp(t[1])] for t in self.int_period])
        year = np.array([[t[0].year, t[1].year] for t in time_array])
        month = np.array([[t[0].month, t[1].month] for t in time_array])
        day = np.array([[t[0].day, t[1].day] for t in time_array])
        doy = np.array([[t[0].timetuple().tm_yday, t[1].timetuple().tm_yday] for t in time_array])
        dtime = np.array([[(t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.), (t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.)] for t in time_array])
        mlat, mlon = self.Apex.geo2apex(self.site[0], self.site[1], self.site[2])
        mlt = np.array([[self.Apex.mlon2mlt(mlon,t[0]), self.Apex.mlon2mlt(mlon,t[1])] for t in time_array])
        return year, month, day, doy, dtime, mlt






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
        chi2 = np.sum((vlos-np.einsum('...i,i->...',A,V))**2/dvlos**2)/(sum(finite)-3)

    except np.linalg.LinAlgError:
        V = np.full(3,np.nan)
        SigV = np.full((3,3),np.nan)
        chi2 = np.nan

    # # if there are too few points for a valid reconstruction, set output to NAN
    # if sum(finite) < minnumpoints:
    #     V = np.full(3,np.nan)
    #     SigV = np.full((3,3),np.nan)

    return V, SigV, chi2


def lin_interp(x, xp, fp, dfp):
    # Piecewise linear interpolation routine that returns interpolated values and their errors

    # find the indicies of xp that bound each value in x
    # Note: where x is out of range of xp, -1 is used as a place holder
    #   This provides a valid "dummy" index for the array calculations and can be used to identify values to nan in final output
    i = np.array([np.argwhere((xi>=xp[:-1]) & (xi<xp[1:])).flatten()[0] if ((xi>=np.nanmin(xp)) & (xi<np.nanmax(xp))) else -1 for xi in x])
    # calculate X
    X = (x-xp[i])/(xp[i+1]-xp[i])
    # calculate interpolated values
    f = (1-X)*fp[i] + X*fp[i+1]
    # calculate interpolation error
    df = np.sqrt((1-X)**2*dfp[i]**2 + X**2*dfp[i+1]**2)
    # replace out-of-range values with NaN
    f[i<0] = np.nan
    df[i<0] = np.nan

    return f, df

def ion_upflow(Te,Ti,ne):
    # calculate ion upflow empirically using the Lu/Zou method (not published yet?)
    pass

def magnitude_direction(A,Sig,e):
    # Calculate the magnitude of vector A and the clockwise angle between vectors e and A
    # Also calculates corresponding errors
    # A = vector
    # Sig = covariance matrix for A
    # e = vector to take the direction relative to
    # ep = e x z (vector perpendicular to e and up)
    # This is all done with somewhat obtuse matrix algebra using einsum to prevent nested for loops
    # Input vectors are assumed to have orthogonal components


    AA = np.einsum('...i,...i->...', A, A)                  # dot product of A and A
    ASA = np.einsum('...i,...ij,...j->...', A, Sig, A)      # matrix multipy A*Sig*A
    ee = np.einsum('...i,...i->...', e, e)                  # dot product of e and e
    eA = np.einsum('...i,...i->...', e, A)                  # dot product of e and A

    # calculate magnitude and magnitude error
    magnitude = np.sqrt(AA)
    mag_err = np.sqrt(ASA/AA)

    # find ep, perpendicular to both e and geodetic up
    ep = np.cross(e,np.array([0,0,1]))
    epep = np.einsum('...i,...i->...', ep, ep)
    epA = np.einsum('...i,...i->...', ep, A)

    B = np.einsum('...ij,...i->...ij',ep,eA)-np.einsum('...ij,...i->...ij',e,epA)     # B = ep(e*A)-e(ep*A) = A x (ep x e)
    BSB = np.einsum('...i,...ij,...j->...', B, Sig, B)      # matrix multipy B*Sig*B

    # calculate direction and direction error
    direction = np.arctan2(np.sqrt(ee)*epA,np.sqrt(epep)*eA)
    dir_err = np.sqrt(epep*ee*BSB)/(ee*epA**2-epep*eA**2)

    return magnitude, mag_err, direction*180./np.pi, dir_err*180./np.pi
