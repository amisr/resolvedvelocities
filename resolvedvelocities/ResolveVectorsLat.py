# ResolveVectors.py

import os
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import socket
import getpass
import platform

import apexpy
import tables
import numpy as np
import datetime as dt
from scipy.spatial import Delaunay

from .marp import Marp
from .plot import summary_plots

import resolvedvelocities as rv


class ResolveVectors(object):
    def __init__(self, config):
        self.configfile = config
        self.read_config(self.configfile)

        # check if output path exists, if not try to create it
        # else raise exception
        self.create_path(self.outfilepath)

        # check if plotting directory exists, if not, create it
        # else raise exception
        if self.plotsavedir:
            self.create_path(self.plotsavedir)


    def create_path(self,path):
        path = os.path.abspath(path)
        os.makedirs(path,exist_ok=True)

    def read_config(self, config_file):

        # read config file
        config = configparser.ConfigParser()
        config.read(config_file)

        # Possibly could done better with converters?  This may by python3 specific though.
        self.datafile = config.get('FILEIO', 'DATAFILE')

        self.outfilename = config.get('FILEIO', 'OUTFILENAME')
        self.chirp = config.getfloat('CONFIG', 'CHIRP')
        self.covar = [float(i) for i in config.get('CONFIG', 'COVAR').split(',')]
        self.altlim = [float(i) for i in config.get('CONFIG', 'ALTLIM').split(',')]
        self.nelim = [float(i) for i in config.get('CONFIG', 'NELIM').split(',')]
        self.chi2lim = [float(i) for i in config.get('CONFIG', 'CHI2LIM').split(',')]
        self.goodfitcode = [float(i) for i in config.get('CONFIG', 'GOODFITCODE').split(',')]
        self.binvert = eval(config.get('CONFIG', 'BINVERT'))
        # can probably change this parameter to a list of start, end, step similar to vvelsAlt
        self.outalt = [float(i) for i in config.get('CONFIG', 'OUTALT').split(',')]
        self.marprot = [float(i) for i in config.get('CONFIG', 'MARPROT').split(',')]

        # optional parameters
        self.plotsavedir = config.get('PLOTTING', 'PLOTSAVEDIR') if config.has_option('PLOTTING', 'PLOTSAVEDIR') else None
        self.upB_beamcode = config.getint('CONFIG', 'UPB_BEAMCODE') if config.has_option('CONFIG', 'UPB_BEAMCODE') else None
        self.ionup = config.get('CONFIG', 'IONUP') if config.has_option('CONFIG', 'IONUP') else None
        self.use_beams = [int(i) for i in config.get('CONFIG', 'USE_BEAMS').split(',')] if config.has_option('CONFIG', 'USE_BEAMS') else None
        self.integration_time = config.getfloat('CONFIG', 'INTTIME') if config.has_option('CONFIG', 'INTTIME') else None
        self.outfilepath = config.get('FILEIO', 'OUTFILEPATH') if config.has_option('FILEIO', 'OUTFILEPATH') else '.'


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

            # ion masses
            ion_mass = infile.get_node('/FittedParams/IonMass')[:]
            nions = len(ion_mass)                                   # number of ions
            ion_idx = np.argwhere(ion_mass==16.).flatten()[0]       # find index for O+

            # line of sight velocity and error
            self.vlos = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(self.alt)))
            self.dvlos = infile.get_node('/FittedParams/Errors')[:,bm_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(self.alt)))

            # chi2 and fitcode (for filtering poor quality data)
            self.chi2 = infile.get_node('/FittedParams/FitInfo/chi2')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))
            self.fitcode = infile.get_node('/FittedParams/FitInfo/fitcode')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))

            # density (for filtering and ion upflow correction)
            self.ne = infile.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))

            # temperature (for ion upflow)
            self.Te = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,-1,1].reshape((len(self.time[:,0]),len(self.alt)))
            Ts = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:nions,1]
            frac = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:nions,0]
            self.Ti = np.sum(Ts*frac,axis=-1).reshape((len(self.time[:,0]),len(self.alt)))

            # get up-B beam velocities for ion outflow correction
            if self.upB_beamcode:
                upB_idx = np.argwhere(self.BeamCodes==self.upB_beamcode).flatten()
                if upB_idx:
                    upB_alt = infile.get_node('/Geomag/Altitude')[upB_idx,:].flatten()
                    upB_vlos = infile.get_node('/FittedParams/Fits')[:,upB_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(upB_alt)))
                    upB_dvlos = infile.get_node('/FittedParams/Errors')[:,upB_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(upB_alt)))
                    self.upB = {'alt':upB_alt, 'vlos':upB_vlos, 'dvlos':upB_dvlos}
                else:
                    print('Warning: upB beam %d not found. Will not perform upflow subtraction.' % self.upB_beamcode)
                    self.upB = None
                    self.ionup = False



    def filter_data(self):
        # filter and adjust data so it is appropriate for Bayesian reconstruction

        with np.errstate(invalid='ignore'):

            # add chirp to LoS velocity
            self.vlos = self.vlos + self.chirp

            # discard data with low density
            inds = np.where((self.ne < self.nelim[0]) | (self.ne > self.nelim[1]))
            self.vlos[inds] = np.nan
            self.dvlos[inds] = np.nan

            # discard data outside of altitude range
            inds = np.where((self.alt < self.altlim[0]*1000.) | (self.alt > self.altlim[1]*1000.))
            self.vlos[:,inds] = np.nan
            self.dvlos[:,inds] = np.nan

            # discard data with extremely high or extremely low chi2 values
            inds = np.where((self.chi2 < self.chi2lim[0]) | (self.chi2 > self.chi2lim[1]))
            self.vlos[inds] = np.nan
            self.dvlos[inds] = np.nan

            # discard data with poor fitcode (fitcodes 1-4 denote solution found, anything else should not be used)
            inds = np.where(~np.isin(self.fitcode, self.goodfitcode))
            self.vlos[inds] = np.nan
            self.vlos[inds] = np.nan


    def transform(self):
        # transform k vectors from geodetic to geomagnetic

        # find indices where nans exist in the altitude array and should be inserted in to other coordinate/component arrays
        replace_nans = np.array([r - i for i,r in enumerate(np.argwhere(np.isnan(self.alt)).flatten())], dtype=int)

        glat = self.lat[np.isfinite(self.lat)]
        glon = self.lon[np.isfinite(self.lon)]
        galt = self.alt[np.isfinite(self.alt)] / 1000.

        # intialize apex coordinates
        self.marp = Marp(date=dt.datetime.utcfromtimestamp(self.time[0,0]), lam0=self.marprot[0], phi0=self.marprot[1])

        # find magnetic latitude and longitude
        # mlat, mlon = self.Apex.geo2apex(glat, glon, galt)
        mlat, mlon = self.marp.geo2marp(glat, glon, galt)
        self.mlat = np.insert(mlat,replace_nans,np.nan)
        self.mlon = np.insert(mlon,replace_nans,np.nan)
        # print(self.mlat, self.mlon)

        # Analogous to apex basis vectors in geodetic coordinates [e n u]
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
                vupflow, dvupflow = lin_interp(self.alt, self.upB['alt'], self.upB['vlos'][t], self.upB['dvlos'][t])
            elif self.ionup == 'EMP':
                # use empirical method to find ion upflow
                # NOTE: NOT DEVELOPED YET!!!
                vupflow, dvupflow = ion_upflow(self.Te, self.Ti, self.ne, self.alt)

            print(vupflow)
            # LoS velocity correction to remove ion upflow
            self.vlos[t] = self.vlos[t] + self.A[:,2]/self.D*vupflow
            # corrected error in new LoS velocities
            self.dvlos[t] = np.sqrt(self.dvlos[t]**2 + self.A[:,2]**2/self.D**2*dvupflow**2)



    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in config file by a list of bin verticies in apex magnetic coordinates
        # the center of each bin is defined as the average of the verticies

        self.bin_mlat = []
        self.bin_mlon = []
        self.bin_idx = []
        for vert in self.binvert:
            vert = np.array(vert)
            hull = Delaunay(vert)

            self.bin_mlat.append(np.nanmean(vert[:,0]))
            self.bin_mlon.append(np.nanmean(vert[:,1]))
            self.bin_idx.append(np.argwhere(hull.find_simplex(np.array([self.mlat, self.mlon]).T)>=0).flatten())

        self.bin_mlat = np.array(self.bin_mlat)
        self.bin_mlon = np.array(self.bin_mlon)



    def get_integration_periods(self):

        if not self.integration_time:
            # if no integration time specified, use original time periods of input files
            self.integration_periods = self.time
            self.integration_indices = range(len(self.time))

        else:
            # if an integration time is given, calculate new time periods
            self.integration_periods = []
            self.integration_indices = []

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
                    self.integration_periods.append([start_time, temp_end_time])
                    self.integration_indices.append(np.array(idx))
                    idx = []
                    start_time = None
                    continue

            self.integration_periods = np.array(self.integration_periods)


    def compute_vector_velocity(self):
        # Iterate over all bins at each integration period and use Heinselman
        # and Nicolls Bayesian reconstruction algorithm to get full vectors

        Velocity = []
        VelocityCovariance = []
        ChiSquared = []
        NumPoints = []

        # For each integration period and bin, calculate covarient components
        # of drift velocity (Ve1, Ve2, Ve3) loop over integration periods
        for tidx in self.integration_indices:
            Vel = []
            SigmaV = []
            Chi2 = []
            NumP = []
            # loop over spatial bins
            for bidx in self.bin_idx:

                # pull out the line of slight measurements for the time period
                # and bins
                vlos = self.vlos[tidx,bidx[:,np.newaxis]].flatten()
                dvlos = self.dvlos[tidx,bidx[:,np.newaxis]].flatten()
                # pull out the k vectors for the bins and duplicate so they
                # match the number of time measurements
                if self.integration_time:
                    A = np.repeat(self.A[bidx], len(tidx), axis=0)
                else:
                    # if no post integraiton, k vectors do not need to be
                    # duplicated
                    A = self.A[bidx]

                # use Heinselman and Nicolls Bayesian reconstruction algorithm
                # to get full vectors
                V, SigV, chi2, num_points = vvels(vlos, dvlos, A, self.covar)

                # append vector and coviarience matrix
                Vel.append(V)
                SigmaV.append(SigV)
                Chi2.append(chi2)
                NumP.append(num_points)

            Velocity.append(Vel)
            VelocityCovariance.append(SigmaV)
            ChiSquared.append(Chi2)
            NumPoints.append(NumP)

        self.Velocity = np.array(Velocity)
        self.VelocityCovariance = np.array(VelocityCovariance)
        self.ChiSquared = np.array(ChiSquared)
        self.NumPoints = np.array(NumPoints)

    def compute_apex_velocity(self):

        d1m, d2m, d3m, e1m, e2m, e3m = self.marp.basevectors_marp(self.bin_mlat,self.bin_mlon,300.,coords='marp')
        self.bin_alat, self.bin_alon = self.marp.marp2apex(self.bin_mlat,self.bin_mlon)
        _,_,_,_,_,_, d1a, d2a, d3a, e1a, e2a, e3a = self.marp.basevectors_apex(self.bin_alat,self.bin_alon,300.,coords='apex')

        em = np.array([e1m, e2m, e3m])
        da = np.array([d1a, d2a, d3a])

        # form transformation array
        R = np.einsum('ik...,jk...->...ij', da, em)

        # transform velocity
        self.Velocity = np.einsum('ijk,...ik->...ij',R,self.Velocity)
        # transfom velocity covarience matrix
        self.VelocityCovariance = np.einsum('ijk,...ikl,iml->...ijm',R,self.VelocityCovariance,R)

    def compute_electric_field(self):
        # calculate electric field

        # find Be3 value at each output bin location
        # NOTE: Be3 is constant along magnetic field lines, so the altitude chosen here doesn't matter
        # Be3, __, __, __ = self.Apex.bvectors_apex(self.bin_mlat,self.bin_mlon,300.,coords='apex')
        Be3, __, __, __ = self.marp.bvectors_apex(self.bin_alat,self.bin_alon,300.,coords='apex')
        # Be3 = np.full(len(self.bin_mlat),1.0)        # set Be3 array to 1.0 - useful for debugging linear algebra

        # form rotation array
        R = np.einsum('i,jk->ijk',Be3,np.array([[0,-1,0],[1,0,0],[0,0,0]]))
        # Calculate contravarient components of electric field (Ed1, Ed2, Ed3)
        self.ElectricField = np.einsum('ijk,...ik->...ij',R,self.Velocity)
        # Calculate electric field covariance matrix (SigE = R*SigV*R.T)
        self.ElectricFieldCovariance = np.einsum('ijk,...ikl,iml->...ijm',R,self.VelocityCovariance,R)


    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        # altitudes are defined by config file

        hbins = len(self.bin_alat)
        vbins = len(self.outalt)
        alat = np.tile(self.bin_alat,vbins)
        alon = np.tile(self.bin_alon,vbins)
        alt = np.repeat(self.outalt, hbins)

        # calculate bin locations in geodetic coordinates
        glat, glon, err = self.marp.apex2geo(alat, alon, alt)
        self.bin_glat = glat.reshape((vbins,hbins))
        self.bin_glon = glon.reshape((vbins,hbins))
        self.bin_galt = alt.reshape((vbins,hbins))

        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.marp.basevectors_apex(glat, glon, alt)

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
        output = os.path.join(self.outfilepath, self.outfilename)
        FILTERS = tables.Filters(complib='zlib', complevel=1)
        with tables.open_file(output, mode='w',filters=FILTERS) as outfile:

            # copy some groups directly from fitted input file
            with tables.open_file(self.datafile, mode='r') as infile:
                outfile.copy_children(infile.get_node('/Site'), outfile.create_group('/','Site'))
                if not self.integration_time:
                    outfile.copy_children(infile.get_node('/Time'), outfile.create_group('/','Time'))
                else:
                    outfile.create_group('/','Time')
                    year, month, day, doy, dtime, mlt = self.create_time_arrays()

                    atom = tables.Atom.from_dtype(self.integration_periods.dtype)
                    arr = outfile.create_carray('/Time', 'UnixTime',atom,self.integration_periods.shape)
                    arr[:] = self.integration_periods
                    outfile.set_node_attr('/Time/UnixTime', 'TITLE', 'UnixTime')
                    outfile.set_node_attr('/Time/UnixTime', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/UnixTime', 'Unit', 'Seconds')

                    atom = tables.Atom.from_dtype(year.dtype)
                    arr = outfile.create_carray('/Time', 'Year',atom,year.shape)
                    arr[:] = year
                    outfile.set_node_attr('/Time/Year', 'TITLE', 'Year')
                    outfile.set_node_attr('/Time/Year', 'Size', 'Nrecords x 2 (Start and end of integration')

                    atom = tables.Atom.from_dtype(month.dtype)
                    arr = outfile.create_carray('/Time', 'Month',atom,month.shape)
                    arr[:] = month
                    outfile.set_node_attr('/Time/Month', 'TITLE', 'Month')
                    outfile.set_node_attr('/Time/Month', 'Size', 'Nrecords x 2 (Start and end of integration')

                    atom = tables.Atom.from_dtype(day.dtype)
                    arr = outfile.create_carray('/Time', 'Day',atom,day.shape)
                    arr[:] = day
                    outfile.set_node_attr('/Time/Day', 'TITLE', 'Day of Month')
                    outfile.set_node_attr('/Time/Day', 'Size', 'Nrecords x 2 (Start and end of integration')

                    atom = tables.Atom.from_dtype(doy.dtype)
                    arr = outfile.create_carray('/Time', 'doy',atom,doy.shape)
                    arr[:] = doy
                    outfile.set_node_attr('/Time/doy', 'TITLE', 'Day of Year')
                    outfile.set_node_attr('/Time/doy', 'Size', 'Nrecords x 2 (Start and end of integration')

                    atom = tables.Atom.from_dtype(dtime.dtype)
                    arr = outfile.create_carray('/Time', 'dtime',atom,dtime.shape)
                    arr[:] = dtime
                    outfile.set_node_attr('/Time/dtime', 'TITLE', 'Decimal Time')
                    outfile.set_node_attr('/Time/dtime', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/dtime', 'Unit', 'UT Hour')

                    atom = tables.Atom.from_dtype(mlt.dtype)
                    arr = outfile.create_carray('/Time', 'MagneticLocalTimeSite',atom,mlt.shape)
                    arr[:] = mlt
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'TITLE', 'Magnetic Local Time')
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'Size', 'Nrecords x 2 (Start and end of integration')
                    outfile.set_node_attr('/Time/MagneticLocalTimeSite', 'Unit', 'UT Hour')

            outfile.create_group('/', 'Magnetic')

            atom = tables.Atom.from_dtype(self.bin_alat.dtype)
            arr = outfile.create_carray('/Magnetic', 'Latitude',atom,self.bin_alat.shape)
            arr[:] = self.bin_alat
            outfile.set_node_attr('/Magnetic/Latitude', 'TITLE', 'Magnetic Latitude')
            outfile.set_node_attr('/Magnetic/Latitude', 'Size', 'Nbins')

            atom = tables.Atom.from_dtype(self.bin_alon.dtype)
            arr = outfile.create_carray('/Magnetic','Longitude',atom,self.bin_alon.shape)
            arr[:] = self.bin_alon
            outfile.set_node_attr('/Magnetic/Longitude', 'TITLE', 'Magnetic Longitude')
            outfile.set_node_attr('/Magnetic/Longitude', 'Size', 'Nbins')

            atom = tables.Atom.from_dtype(self.Velocity.dtype)
            arr = outfile.create_carray('/Magnetic', 'Velocity',atom,self.Velocity.shape)
            arr[:] = self.Velocity
            outfile.set_node_attr('/Magnetic/Velocity', 'TITLE', 'Plasma Drift Velocity')
            outfile.set_node_attr('/Magnetic/Velocity', 'Size', 'Nrecords x Nbins x 3 (Ve1, Ve2, Ve3)')
            outfile.set_node_attr('/Magnetic/Velocity', 'Units', 'm/s')

            atom = tables.Atom.from_dtype(self.VelocityCovariance.dtype)
            arr = outfile.create_carray('/Magnetic','CovarianceV',atom,self.VelocityCovariance.shape)
            arr[:] = self.VelocityCovariance
            outfile.set_node_attr('/Magnetic/CovarianceV', 'TITLE', 'Velocity Covariance Matrix')
            outfile.set_node_attr('/Magnetic/CovarianceV', 'Size', 'Nrecords x Nbins x 3 x 3')
            outfile.set_node_attr('/Magnetic/CovarianceV', 'Units', '(m/s)^2')

            atom = tables.Atom.from_dtype(self.ElectricField.dtype)
            arr = outfile.create_carray('/Magnetic','ElectricField',atom,self.ElectricField.shape)
            arr[:] = self.ElectricField
            outfile.set_node_attr('/Magnetic/ElectricField', 'TITLE', 'Convection Electric Field')
            outfile.set_node_attr('/Magnetic/ElectricField', 'Size', 'Nrecords x Nbins x 3 (Ed1, Ed2, Ed3)')
            outfile.set_node_attr('/Magnetic/ElectricField', 'Units', 'V/m')

            atom = tables.Atom.from_dtype(self.ElectricFieldCovariance.dtype)
            arr = outfile.create_carray('/Magnetic','CovarianceE',atom,self.ElectricFieldCovariance.shape)
            arr[:] = self.ElectricFieldCovariance
            outfile.set_node_attr('/Magnetic/CovarianceE', 'TITLE', 'Electric Field Covariance Matrix')
            outfile.set_node_attr('/Magnetic/CovarianceE', 'Size', 'Nrecords x Nbins x 3 x 3')
            outfile.set_node_attr('/Magnetic/CovarianceE', 'Units', '(V/m)^2')


            outfile.create_group('/', 'Geodetic')

            atom = tables.Atom.from_dtype(self.bin_glat.dtype)
            arr = outfile.create_carray('/Geodetic', 'Latitude',atom,self.bin_glat.shape)
            arr[:] = self.bin_glat
            outfile.set_node_attr('/Geodetic/Latitude', 'TITLE', 'Geographic Latitude')
            outfile.set_node_attr('/Geodetic/Latitude', 'Size', 'Nalt x Nbins')

            atom = tables.Atom.from_dtype(self.bin_glon.dtype)
            arr = outfile.create_carray('/Geodetic','Longitude',atom,self.bin_glon.shape)
            arr[:] = self.bin_glon
            outfile.set_node_attr('/Geodetic/Longitude', 'TITLE', 'Geographic Longitude')
            outfile.set_node_attr('/Geodetic/Longitude', 'Size', 'Nalt x Nbins')

            atom = tables.Atom.from_dtype(self.bin_galt.dtype)
            arr = outfile.create_carray('/Geodetic','Altitude',atom,self.bin_galt.shape)
            arr[:] = self.bin_galt
            outfile.set_node_attr('/Geodetic/Altitude', 'TITLE', 'Geographic Altitude')
            outfile.set_node_attr('/Geodetic/Altitude', 'Size', 'Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/Altitude', 'Units', 'km')

            atom = tables.Atom.from_dtype(self.Velocity_gd.dtype)
            arr = outfile.create_carray('/Geodetic', 'Velocity',atom,self.Velocity_gd.shape)
            arr[:] = self.Velocity_gd
            outfile.set_node_attr('/Geodetic/Velocity', 'TITLE', 'Plasma Drift Velocity')
            outfile.set_node_attr('/Geodetic/Velocity', 'Size', 'Nrecords x Nalt x Nbins x 3 (East, North, Up)')
            outfile.set_node_attr('/Geodetic/Velocity', 'Units', 'm/s')

            atom = tables.Atom.from_dtype(self.VelocityCovariance_gd.dtype)
            arr = outfile.create_carray('/Geodetic','CovarianceV',atom,self.VelocityCovariance_gd.shape)
            arr[:] = self.VelocityCovariance_gd
            outfile.set_node_attr('/Geodetic/CovarianceV', 'TITLE', 'Velocity Covariance Matrix')
            outfile.set_node_attr('/Geodetic/CovarianceV', 'Size', 'Nrecords x Nalt x Nbins x 3 x 3')
            outfile.set_node_attr('/Geodetic/CovarianceV', 'Units', '(m/s)^2')

            atom = tables.Atom.from_dtype(self.Vgd_mag.dtype)
            arr = outfile.create_carray('/Geodetic','Vmag',atom,self.Vgd_mag.shape)
            arr[:] = self.Vgd_mag
            outfile.set_node_attr('/Geodetic/Vmag', 'TITLE', 'Velocity Magnitude')
            outfile.set_node_attr('/Geodetic/Vmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/Vmag', 'Units', 'm/s')

            atom = tables.Atom.from_dtype(self.Vgd_mag_err.dtype)
            arr = outfile.create_carray('/Geodetic','errVmag',atom,self.Vgd_mag_err.shape)
            arr[:] = self.Vgd_mag_err
            outfile.set_node_attr('/Geodetic/errVmag', 'TITLE', 'Velocity Magnitude Error')
            outfile.set_node_attr('/Geodetic/errVmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/errVmag', 'Units', 'm/s')

            atom = tables.Atom.from_dtype(self.Vgd_dir.dtype)
            arr = outfile.create_carray('/Geodetic','Vdir',atom,self.Vgd_dir.shape)
            arr[:] = self.Vgd_dir
            outfile.set_node_attr('/Geodetic/Vdir', 'TITLE', 'Velocity Direction Angle East of North Magnetic Meridian (-e2)')
            outfile.set_node_attr('/Geodetic/Vdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/Vdir', 'Units', 'Degrees')

            atom = tables.Atom.from_dtype(self.Vgd_dir_err.dtype)
            arr = outfile.create_carray('/Geodetic','errVdir',atom,self.Vgd_dir_err.shape)
            arr[:] = self.Vgd_dir_err
            outfile.set_node_attr('/Geodetic/errVdir', 'TITLE', 'Error in Velocity Direction')
            outfile.set_node_attr('/Geodetic/errVdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/errVdir', 'Units', 'Degrees')

            atom = tables.Atom.from_dtype(self.ElectricField_gd.dtype)
            arr = outfile.create_carray('/Geodetic','ElectricField',atom,self.ElectricField_gd.shape)
            arr[:] = self.ElectricField_gd
            outfile.set_node_attr('/Geodetic/ElectricField', 'TITLE', 'Convection Electric Field')
            outfile.set_node_attr('/Geodetic/ElectricField', 'Size', 'Nrecords x Nalt x Nbins x 3 (East, North, Up)')
            outfile.set_node_attr('/Geodetic/ElectricField', 'Units', 'V/m')

            atom = tables.Atom.from_dtype(self.ElectricFieldCovariance_gd.dtype)
            arr = outfile.create_carray('/Geodetic','CovarianceE',atom,self.ElectricFieldCovariance_gd.shape)
            arr[:] = self.ElectricFieldCovariance_gd
            outfile.set_node_attr('/Geodetic/CovarianceE', 'TITLE', 'Electric Field Covariance Matrix')
            outfile.set_node_attr('/Geodetic/CovarianceE', 'Size', 'Nrecords x Nalt x Nbins x 3 x 3')
            outfile.set_node_attr('/Geodetic/CovarianceE', 'Units', '(V/m)^2')

            atom = tables.Atom.from_dtype(self.Egd_mag.dtype)
            arr = outfile.create_carray('/Geodetic','Emag',atom,self.Egd_mag.shape)
            arr[:] = self.Egd_mag
            outfile.set_node_attr('/Geodetic/Emag', 'TITLE', 'Electric Field Magnitude')
            outfile.set_node_attr('/Geodetic/Emag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/Emag', 'Units', 'V/m')

            atom = tables.Atom.from_dtype(self.Egd_mag_err.dtype)
            arr = outfile.create_carray('/Geodetic','errEmag',atom,self.Egd_mag_err.shape)
            arr[:] = self.Egd_mag_err
            outfile.set_node_attr('/Geodetic/errEmag', 'TITLE', 'Electric Field Magnitude Error')
            outfile.set_node_attr('/Geodetic/errEmag', 'Size', 'Nrecords x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/errEmag', 'Units', 'V/m')

            atom = tables.Atom.from_dtype(self.Egd_dir.dtype)
            arr = outfile.create_carray('/Geodetic','Edir',atom,self.Egd_dir.shape)
            arr[:] = self.Egd_dir
            outfile.set_node_attr('/Geodetic/Edir', 'TITLE', 'Electric Field Direction Angle East of North Magnetic Meridian (-e2)')
            outfile.set_node_attr('/Geodetic/Edir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/Edir', 'Units', 'Degrees')

            atom = tables.Atom.from_dtype(self.Egd_dir_err.dtype)
            arr = outfile.create_carray('/Geodetic','errEdir',atom,self.Egd_dir_err.shape)
            arr[:] = self.Egd_dir_err
            outfile.set_node_attr('/Geodetic/errEdir', 'TITLE', 'Error in Electric Field Direction')
            outfile.set_node_attr('/Geodetic/errEdir', 'Size', 'Nrecord x Nalt x Nbins')
            outfile.set_node_attr('/Geodetic/errEdir', 'Units', 'Degrees')

            outfile.create_group('/', 'ProcessingParams')

            atom = tables.Atom.from_dtype(self.ChiSquared.dtype)
            arr = outfile.create_carray('/ProcessingParams', 'Chi2',atom,self.ChiSquared.shape)
            arr[:] = self.ChiSquared
            outfile.set_node_attr('/ProcessingParams/Chi2', 'TITLE', 'Reduced Chi-Squared')
            outfile.set_node_attr('/ProcessingParams/Chi2', 'Size', 'Nrecords x Nbins')

            atom = tables.Atom.from_dtype(self.NumPoints.dtype)
            arr = outfile.create_carray('/ProcessingParams', 'NumPoints',atom,self.NumPoints.shape)
            arr[:] = self.NumPoints
            outfile.set_node_attr('/ProcessingParams/NumPoints', 'TITLE', 'Number of input data points used to estimate the vector')
            outfile.set_node_attr('/ProcessingParams/NumPoints', 'Size', 'Nrecords x Nbins')

            outfile.create_group('/ProcessingParams', 'Apexpy')

            outfile.create_array('/ProcessingParams/Apexpy', 'Year', self.marp.year)
            outfile.set_node_attr('/ProcessingParams/Apexpy/Year', 'TITLE', 'Decimal Year used for IGRF Model')

            outfile.create_array('/ProcessingParams/Apexpy','RefHeight',self.marp.refh)
            outfile.set_node_attr('/ProcessingParams/Apexpy/RefHeight', 'TITLE', 'Reference height used for Apex coordinates')
            outfile.set_node_attr('/ProcessingParams/Apexpy/RefHeight', 'Units', 'km')

            outfile.create_array('/ProcessingParams/Apexpy','Version',apexpy.__version__.encode('utf-8'))
            outfile.set_node_attr('/ProcessingParams/Apexpy/Version', 'TITLE', 'Apexpy version used.')

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

            # config file
            outfile.create_group('/ProcessingParams', 'ConfigFile')

            outfile.create_array('/ProcessingParams/ConfigFile', 'Name', Name.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ConfigFile', 'Path', Path.encode('utf-8'))
            outfile.create_array('/ProcessingParams/ConfigFile', 'Contents', Contents.encode('utf-8'))

            # resolved velocities code version
            outfile.create_array('/ProcessingParams', 'SoftwareVersion', rv.__version__.encode('utf-8'))
            outfile.set_node_attr('/ProcessingParams/SoftwareVersion', 'TITLE', 'Version of the resolvedvelocities software that was used.')

            # input file
            outfile.create_array('/ProcessingParams', 'InputFile', self.datafile.encode('utf-8'))




    def create_time_arrays(self):
        time_array = np.array([[dt.datetime.utcfromtimestamp(t[0]), dt.datetime.utcfromtimestamp(t[1])] for t in self.integration_periods])
        year = np.array([[t[0].year, t[1].year] for t in time_array])
        month = np.array([[t[0].month, t[1].month] for t in time_array])
        day = np.array([[t[0].day, t[1].day] for t in time_array])
        doy = np.array([[t[0].timetuple().tm_yday, t[1].timetuple().tm_yday] for t in time_array])
        dtime = np.array([[(t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.), (t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.)] for t in time_array])
        mlat, mlon = self.Apex.geo2apex(self.site[0], self.site[1], self.site[2])
        mlt = np.array([[self.Apex.mlon2mlt(mlon,t[0]), self.Apex.mlon2mlt(mlon,t[1])] for t in time_array])
        return year, month, day, doy, dtime, mlt


    def create_plots(self,alt=300.,vcomptitles=None,vcomplim=None,
                     vcompcmap=None,ecomptitles=None,ecomplim=None,
                     ecompcmap=None,vmagtitles=None,vmaglim=None,
                     vmagcmap=None,emagtitles=None,emaglim=None,
                     emagcmap=None):

        if self.plotsavedir:

            # break up arrays into chunks of time no bigger than 24 hours
            chunks_to_plot = list()

            start_time = None
            start_ind = None
            num_times = len(self.integration_periods)
            for i,time_pair in enumerate(self.integration_periods):
                temp_start_time, temp_end_time = time_pair
                if start_time is None:
                    start_time = temp_start_time
                    start_ind = i
                time_diff = temp_end_time - start_time

                if (time_diff >= 24*3600) or (i == num_times -1):
                    chunks_to_plot.append([start_ind,i])
                    start_time = None
                    start_ind = None
                    continue

            num_chunks = len(chunks_to_plot)
            for t, [start_ind,end_ind] in enumerate(chunks_to_plot):
                # if only 1 day worth of data, set t=None so we don't have a
                # 'byDay' in the plot file names
                if (num_chunks == 1):
                    t = None

                # make vector plots
                times = self.integration_periods[start_ind:end_ind,:]
                vels = self.Velocity[start_ind:end_ind,:]
                covvels = self.VelocityCovariance[start_ind:end_ind,:]
                summary_plots.plot_components(times,self.bin_mlat,
                                              self.bin_mlon,vels,covvels,
                                              param='V',titles=vcomptitles,
                                              clim=vcomplim,cmap=vcompcmap,
                                              savedir=self.plotsavedir,
                                              savenamebase=self.outfilename,
                                              day_ind=t)
                es = self.ElectricField[start_ind:end_ind,:]
                coves = self.ElectricFieldCovariance[start_ind:end_ind,:]
                summary_plots.plot_components(times,self.bin_mlat,
                                              self.bin_mlon,es,coves,
                                              param='E',titles=ecomptitles,
                                              clim=ecomplim,cmap=ecompcmap,
                                              savedir=self.plotsavedir,
                                              savenamebase=self.outfilename,
                                              day_ind=t)

                # make magnitude plots
                # find index of altitude bin that is closest to alt
                i = np.argmin(np.abs(self.bin_galt[:,0]-alt))
                vmag = self.Vgd_mag[start_ind:end_ind,i,:]
                evmag = self.Vgd_mag_err[start_ind:end_ind,i,:]
                vdir = self.Vgd_dir[start_ind:end_ind,i,:]
                evdir = self.Vgd_dir_err[start_ind:end_ind,i,:]
                chi2 = self.ChiSquared[start_ind:end_ind,:]
                summary_plots.plot_magnitude(times,self.bin_mlat,self.bin_mlon,
                                             vmag,evmag,vdir,evdir,chi2,
                                             param='V',titles=vmagtitles,
                                             clim=vmaglim,cmap=vmagcmap,
                                             savedir=self.plotsavedir,
                                             savenamebase=self.outfilename,
                                             day_ind=t)

                emag = self.Egd_mag[start_ind:end_ind,i,:]
                eemag = self.Egd_mag_err[start_ind:end_ind,i,:]
                edir = self.Egd_dir[start_ind:end_ind,i,:]
                eedir = self.Egd_dir_err[start_ind:end_ind,i,:]
                summary_plots.plot_magnitude(times,self.bin_mlat,self.bin_mlon,
                                             emag,eemag,edir,eedir,chi2,
                                             param='E',titles=emagtitles,
                                             clim=emaglim,cmap=emagcmap,
                                             savedir=self.plotsavedir,
                                             savenamebase=self.outfilename,
                                             day_ind=t)




def vvels(vlos, dvlos, A, cov, minnumpoints=1):
    # Bayesian inference method described in Heinselman and Nicolls 2008
    # vector velocity algorithm

    # Get indices for finite valued data points
    finite = np.isfinite(vlos)
    num_points = np.sum(finite)
    dof = num_points - 3 # solving for 3 components

    # Filter inputs to only use finite valued data
    vlos = vlos[finite]
    dvlos = dvlos[finite]
    A = A[finite]

    SigmaE = np.diagflat(dvlos**2)
    SigmaV = np.diagflat(cov)

    try:
        # measurement errors and a priori covariance, terms in the inverse
        # (Heinselman and Nicolls 2008 eqn 12)
        # I = (A*SigV*A.T + SigE)^-1
        I = np.linalg.inv(np.einsum('jk,kl,ml->jm',A,SigmaV,A) + SigmaE)
        # calculate velocity estimate (Heinselman and Nicolls 2008 eqn 12)
        V = np.einsum('jk,lk,lm,m->j',SigmaV,A,I,vlos)
        # calculate covariance of velocity estimate
        # (Heinselman and Nicolls 2008 eqn 13)
        term1 = np.einsum('kj,kl,lm->jm',A,np.linalg.inv(SigmaE),A)
        term2 = np.linalg.inv(SigmaV)
        SigV = np.linalg.inv(term1 + term2)

        # chi-squared is meaningless for an underdetermined problem
        if dof < 1:
            chi2 = np.nan
        else:
            model = np.einsum('...i,i->...',A,V)
            chi2 = np.sum((vlos - model)**2 / dvlos**2) / dof

    except np.linalg.LinAlgError:
        V = np.full(3,np.nan)
        SigV = np.full((3,3),np.nan)
        chi2 = np.nan
        num_points = 0

    return V, SigV, chi2, num_points


def lin_interp(x, xp, fp, dfp):
    # Piecewise linear interpolation routine that returns interpolated values
    # and their errors

    # find the indicies of xp that bound each value in x
    # Note: where x is out of range of xp, -1 is used as a place holder
    #   This provides a valid "dummy" index for the array calculations and can
    # be used to identify values to nan in final output
    xpmin = np.nanmin(xp)
    xpmax = np.nanmax(xp)
    i = np.array([np.argwhere((xi>=xp[:-1]) & (xi<xp[1:])).flatten()[0]
                  if ((xi>=xpmin) & (xi<xpmax)) else -1 for xi in x])
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
    # Calculate the magnitude of vector A and the clockwise angle between
    # vectors e and A
    # Also calculates corresponding errors
    # A = vector
    # Sig = covariance matrix for A
    # e = vector to take the direction relative to
    # ep = e x z (vector perpendicular to e and up)
    # This is all done with an error analysis using addition in quadrature.
    # All calculations are done with matrix algebra using einsum to prevent
    # nested for loops.
    # Input vectors are assumed to have orthogonal components

    # Calculate the Magnitude of input A

    # Some helper matricies
    # dot product of A and A
    AA = np.einsum('...i,...i->...', A, A)
    # matrix multipy A*Sig*transpose(A)
    ASA = np.einsum('...i,...ij,...j->...', A, Sig, A)

    # calculate magnitude and magnitude error
    magnitude = np.sqrt(AA)
    mag_err = np.sqrt(ASA/AA)


    # Now find the angle clockwise around geodetic up in the horizontal plane
    # that is perpendicular to geodetic up, where angle=0 is along the
    # projection 'e' in the horizontal plane

    # Some helper matricies
    # dot product of e and e
    ee = np.einsum('...i,...i->...', e, e)
    # dot product of e and A
    eA = np.einsum('...i,...i->...', e, A)
    # find ep, perpendicular to both e and geodetic up
    ep = np.cross(e,np.array([0,0,1]))
    epep = np.einsum('...i,...i->...', ep, ep)
    epA = np.einsum('...i,...i->...', ep, A)

    # B = ep(e*A)-e(ep*A) = A x (ep x e)
    B = np.einsum('...ij,...i->...ij',ep,eA)-np.einsum('...ij,...i->...ij',e,epA)
    # matrix multipy B*Sig*B (covariance propagation)
    BSB = np.einsum('...i,...ij,...j->...', B, Sig, B)

    # calculate direction and direction error
    direction = np.arctan2(np.sqrt(ee)*epA,np.sqrt(epep)*eA)
    dir_err = np.sqrt(epep*ee*BSB)/(ee*epA**2-epep*eA**2)

    return magnitude, mag_err, direction*180./np.pi, dir_err*180./np.pi
