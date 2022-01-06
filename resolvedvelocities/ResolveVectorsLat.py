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

from .DataHandler import FittedVelocityDataHandler
from .marp import Marp
from .plot import summary_plots
from .utils import *
import resolvedvelocities as rv


class ResolveVectorsLat(object):
    def __init__(self, configfile):
        self.configfile = configfile
        self.read_config(self.configfile)

        print(self.datafile)
        self.datahandler = FittedVelocityDataHandler(self.datafile)
        self.datahandler.load_data(self.use_beams)
        self.datahandler.filter(chi2=self.chi2lim, ne=self.nelim, alt=self.altlim, fitcode=self.goodfitcode, chirp=self.chirp)

        if self.integration_time:
            self.integration_periods = get_integration_periods(self.datahandler.utime, self.integration_time)
        else:
            self.integration_periods = self.datahandler.utime

        self.outalt = self.create_alt_array()



    def run(self):

        self.transform()
        self.bin_data()
        self.compute_vector_velocity()
        self.compute_apex_velocity()
        self.compute_electric_field()
        self.compute_geodetic_output()
        self.save_output()
        if self.plotsavedir:
            self.create_plots()



    def read_config(self, config_file):

        # read config file
        config = configparser.ConfigParser(converters={'list':parse_list})
        config.read(config_file)

        self.datafile = config.get('FILEIO', 'DATAFILE')

        self.outfilename = config.get('FILEIO', 'OUTFILENAME')
        self.chirp = config.getfloat('CONFIG', 'CHIRP')
        self.covar = config.getlist('CONFIG', 'COVAR')
        self.altlim = config.getlist('CONFIG', 'ALTLIM')
        self.nelim = config.getlist('CONFIG', 'NELIM')
        self.chi2lim = config.getlist('CONFIG', 'CHI2LIM')
        self.goodfitcode = config.getlist('CONFIG', 'GOODFITCODE')
        self.binvert = eval(config.get('CONFIG', 'BINVERT'))
        self.altitude_bins_def = config.get('CONFIG', 'OUTALT')

        # optional parameters
        self.marprot = config.getlist('CONFIG', 'MARPROT') if config.has_option('CONFIG', 'MARPROT') else [0.,0.]
        self.plotsavedir = config.get('PLOTTING', 'PLOTSAVEDIR') if config.has_option('PLOTTING', 'PLOTSAVEDIR') else None
        self.upB_beamcode = config.getint('CONFIG', 'UPB_BEAMCODE') if config.has_option('CONFIG', 'UPB_BEAMCODE') else None
        self.ionup = config.get('CONFIG', 'IONUP') if config.has_option('CONFIG', 'IONUP') else None
        self.use_beams = config.getlist('CONFIG', 'USE_BEAMS') if config.has_option('CONFIG', 'USE_BEAMS') else None
        self.integration_time = config.getfloat('CONFIG', 'INTTIME') if config.has_option('CONFIG', 'INTTIME') else None
        self.outfilepath = config.get('FILEIO', 'OUTFILEPATH') if config.has_option('FILEIO', 'OUTFILEPATH') else '.'


    def create_alt_array(self):

        altarr = np.empty((0,))
        groups = self.altitude_bins_def.split(';')
        for i,group in enumerate(groups):
            start, stop, step = [float(i) for i in group.split(',')]
            altarr = np.append(altarr, np.arange(start, stop, step))
        return altarr


    def transform(self):
        # transform k vectors from geodetic to geomagnetic

        # intialize apex coordinates
        self.marp = Marp(date=dt.datetime.utcfromtimestamp(self.datahandler.utime[0,0]), lam0=self.marprot[0], phi0=self.marprot[1])

        # find magnetic latitude and longitude
        self.mlat, self.mlon = self.marp.geo2marp(self.datahandler.lat, self.datahandler.lon, self.datahandler.alt)

        # Analogous to apex basis vectors in geodetic coordinates [e n u]
        # f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, galt)
        d1, d2, d3, e1, e2, e3 = self.marp.basevectors_marp(self.datahandler.lat, self.datahandler.lon, self.datahandler.alt)
        e = np.array([e1,e2,e3]).T

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([self.datahandler.ke, self.datahandler.kn, self.datahandler.kz]).T

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
                vupflow, dvupflow = lin_interp(self.datahandler.alt, self.upB['alt'], self.upB['vlos'][t], self.upB['dvlos'][t])
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




    def compute_vector_velocity(self):
        # Iterate over all bins at each integration period and use Heinselman
        # and Nicolls Bayesian reconstruction algorithm to get full vectors

        num_integrations = len(self.integration_periods)
        num_bins = len(self.bin_mlat)

        self.Velocity = np.full((num_integrations,num_bins,3), np.nan)
        self.VelocityCovariance = np.full((num_integrations,num_bins,3,3), np.nan)
        self.ChiSquared = np.full((num_integrations,num_bins), np.nan)
        self.NumPoints = np.full((num_integrations,num_bins), np.nan)

        # For each integration period and bin, calculate covarient components
        # of drift velocity (Ve1, Ve2, Ve3) loop over integration periods
        for i, integration_period in enumerate(self.integration_periods):

            tidx = self.datahandler.get_record_indices(integration_period[0],integration_period[1])

            # loop over spatial bins
            for k, bidx in enumerate(self.bin_idx):

                # pull out the line of slight measurements for the time period
                # and bins
                vlos = self.datahandler.vlos[tidx,bidx[:,np.newaxis]].flatten()
                dvlos = self.datahandler.dvlos[tidx,bidx[:,np.newaxis]].flatten()
                A = np.repeat(self.A[bidx], len(tidx), axis=0)

                # use Heinselman and Nicolls Bayesian reconstruction algorithm
                # to get full vectors
                V, SigV, chi2, N = vvels(vlos, dvlos, A, self.covar)

                self.Velocity[i,k,:] = V
                self.VelocityCovariance[i,k,:] = SigV
                self.ChiSquared[i,k] = chi2
                self.NumPoints[i,k] = N



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
        # save output file

        os.makedirs(os.path.abspath(self.outfilepath),exist_ok=True)
        output = os.path.join(self.outfilepath, self.outfilename)
        FILTERS = tables.Filters(complib='zlib', complevel=1)
        with tables.open_file(output, mode='w',filters=FILTERS) as outfile:

            # copy some groups directly from fitted input file
            with tables.open_file(self.datafile, mode='r') as infile:
                outfile.copy_children(infile.get_node('/Site'), outfile.create_group('/','Site'))

            outfile.create_group('/','Time')
            year, month, day, doy, dtime, mlt = create_time_arrays(self.integration_periods, self.datahandler.site)

            save_carray(outfile, '/Time/UnixTime', self.integration_periods, {'TITLE':'UnixTime', 'Size':'Nrecords x 2 (Start and end of integration)', 'Units':'Seconds'})
            save_carray(outfile, '/Time/Year', year, {'TITLE':'Year', 'Size':'Nrecords x 2 (Start and end of integration)'})
            save_carray(outfile, '/Time/Month', month, {'TITLE':'Month', 'Size':'Nrecords x 2 (Start and end of integration)'})
            save_carray(outfile, '/Time/Day', day, {'TITLE':'Day of Month', 'Size':'Nrecords x 2 (Start and end of integration)'})
            save_carray(outfile, '/Time/doy', doy, {'TITLE':'Day of Year', 'Size':'Nrecords x 2 (Start and end of integration)'})
            save_carray(outfile, '/Time/dtime', dtime, {'TITLE':'Decimal Hour of Day', 'Size':'Nrecords x 2 (Start and end of integration)'})
            save_carray(outfile, '/Time/MagneticLocalTimeSite', mlt, {'TITLE':'Magnetic Local Time of Site', 'Size':'Nrecords x 2 (Start and end of integration)'})


            outfile.create_group('/', 'Magnetic')

            save_carray(outfile, '/Magnetic/Latitude', self.bin_alat, {'TITLE':'Magnetic Latitude', 'Size':'Nbins'})
            save_carray(outfile, '/Magnetic/Longitude', self.bin_alon, {'TITLE':'Magnetic Longitude', 'Size':'Nbins'})
            save_carray(outfile, '/Magnetic/Velocity', self.Velocity, {'TITLE':'Plasma Drift Velocity', 'Size':'Nrecords x Nbins x 3 (Ve1, Ve2, Ve3)', 'Units':'m/s'})
            save_carray(outfile, '/Magnetic/CovarianceV', self.VelocityCovariance, {'TITLE':'Velocity Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(m/s)^2'})
            save_carray(outfile, '/Magnetic/ElectricField', self.ElectricField, {'TITLE':'Convection Electric Field', 'Size':'Nrecords x Nbins x 3 (Ed1, Ed2, Ed3)', 'Units':'V/m'})
            save_carray(outfile, '/Magnetic/CovarianceE', self.ElectricFieldCovariance, {'TITLE':'Electric Field Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(V/m)^2'})



            outfile.create_group('/', 'Geodetic')

            save_carray(outfile, '/Geodetic/Latitude', self.bin_alat, {'TITLE':'Geographic Latitude', 'Size':'Nalts x Nbins'})
            save_carray(outfile, '/Geodetic/Longitude', self.bin_alon, {'TITLE':'Geographic Longitude', 'Size':'Nalts x Nbins'})
            save_carray(outfile, '/Geodetic/Altitude', self.outalt, {'TITLE':'Geographic Altitude', 'Size':'Nalts x Nbins', 'Units':'km'})
            save_carray(outfile, '/Geodetic/Velocity', self.Velocity_gd, {'TITLE':'Plasma Drift Velocity', 'Size':'Nrecords x Nalts x Nbins x 3 (East, North, Up)', 'Units':'m/s'})
            save_carray(outfile, '/Geodetic/CovarianceV', self.VelocityCovariance_gd, {'TITLE':'Velocity Covariance Matrix', 'Size':'Nrecords x Nalts x Nbins x 3 x 3', 'Units':'(m/s)^2'})
            save_carray(outfile, '/Geodetic/Vmag', self.Vgd_mag, {'TITLE':'Velocity Magnitude', 'Size':'Nrecords x Nalts x Nbins', 'Units':'m/s'})
            save_carray(outfile, '/Geodetic/errVmag', self.Vgd_mag_err, {'TITLE':'Velocity Magnitude Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'m/s'})
            save_carray(outfile, '/Geodetic/Vdir', self.Vgd_dir, {'TITLE':'Velocity Direction Angle East of North Magnetic Meridian (-e2)', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            save_carray(outfile, '/Geodetic/errVdir', self.Vgd_dir_err, {'TITLE':'Velocity Direction Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            save_carray(outfile, '/Geodetic/ElectricField', self.ElectricField_gd, {'TITLE':'Convection Electric Field', 'Size':'Nrecords x Nbins x 3 (East, North, Up)', 'Units':'V/m'})
            save_carray(outfile, '/Geodetic/CovarianceE', self.ElectricFieldCovariance_gd, {'TITLE':'Electric Field Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(V/m)^2'})
            save_carray(outfile, '/Geodetic/Emag', self.Egd_mag, {'TITLE':'Electric Field Magnitude', 'Size':'Nrecords x Nalts x Nbins', 'Units':'V/m'})
            save_carray(outfile, '/Geodetic/errEmag', self.Egd_mag_err, {'TITLE':'Electric Field Magnitude Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'V/m'})
            save_carray(outfile, '/Geodetic/Edir', self.Egd_dir, {'TITLE':'Electric Field Direction Angle East of North Magnetic Meridian (-e2)', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            save_carray(outfile, '/Geodetic/errEdir', self.Egd_dir_err, {'TITLE':'Electric Field Direction Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})


            outfile.create_group('/', 'ProcessingParams')

            save_carray(outfile, '/ProcessingParams/Chi2', self.ChiSquared, {'TITLE':'Reduced Chi-Squared', 'Size':'Nrecords x Nbins'})
            save_carray(outfile, '/ProcessingParams/NumPoints', self.NumPoints, {'TITLE':'Number of input data points used to estimate the vector', 'Size':'Nrecords x Nbins'})

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




    def create_plots(self, alt=300.):

        os.makedirs(os.path.abspath(self.plotsavedir),exist_ok=True)

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
                vcom_fname = 'velocity_components.png'
                ecom_fname = 'electric_field_components.png'
                vmag_fname = 'velocity_magnitudes.png'
                emag_fname = 'electric_field_magnitudes.png'
            else:
                vcom_fname = 'velocity_components_{}.png'.format(t)
                ecom_fname = 'electric_field_components_{}.png'.format(t)
                vmag_fname = 'velocity_magnitudes_{}.png'.format(t)
                emag_fname = 'electric_field_magnitudes_{}.png'.format(t)

            # make vector plots
            times = self.integration_periods[start_ind:end_ind,:]
            binmlat = ['{:.2f} N'.format(ml) for ml in self.bin_mlat]

            vels = self.Velocity[start_ind:end_ind,:]
            covvels = self.VelocityCovariance[start_ind:end_ind,:]
            efs = self.ElectricField[start_ind:end_ind,:]*1000.
            covefs = self.ElectricFieldCovariance[start_ind:end_ind,:]*1000.*1000.

            summary_plots.plot_components(times, self.bin_mlat, vels, covvels,
                            titles=['Ve1 (m/s)','Ve2 (m/s)','Ve3 (m/s)'],
                            ylabel='Apex MLAT', yticklabels=binmlat,
                            clim=[[-1500.,1500.], [0.,350.]], cmap=['coolwarm', 'turbo'],
                            filename=os.path.join(self.plotsavedir,vcom_fname), scale_factors=[1,1,10])

            summary_plots.plot_components(times, self.bin_mlat, efs, covefs,
                            titles=['Ed1 (mV/m)','Ed2 (mV/m)','Ed3 (mV/m)'],
                            ylabel='Apex MLAT', yticklabels=binmlat,
                            clim=[[-75., 75.], [0., 15.]], cmap=['coolwarm', 'turbo'],
                            filename=os.path.join(self.plotsavedir,ecom_fname), scale_factors=[1,1,10])



            # make magnitude plots
            # find index of altitude bin that is closest to alt
            i = np.argmin(np.abs(self.bin_galt[:,0]-alt))
            vmag = self.Vgd_mag[start_ind:end_ind,i,:]
            dvmag = self.Vgd_mag_err[start_ind:end_ind,i,:]
            vdir = self.Vgd_dir[start_ind:end_ind,i,:]
            dvdir = self.Vgd_dir_err[start_ind:end_ind,i,:]
            emag = self.Egd_mag[start_ind:end_ind,i,:]*1000.
            demag = self.Egd_mag_err[start_ind:end_ind,i,:]*1000.
            edir = self.Egd_dir[start_ind:end_ind,i,:]
            dedir = self.Egd_dir_err[start_ind:end_ind,i,:]
            chi2 = self.ChiSquared[start_ind:end_ind,:]

            titles = ['V mag. (m/s)', 'V mag. err. (m/s)', 'V dir. (deg)', 'V dir. err. (deg)', '']
            clim = [[0.,1500.],[0., 350.],[-180., 180.],[0., 35.]]
            cmap = ['viridis', 'turbo', 'twilight', 'turbo']

            summary_plots.plot_magnitude(times, self.bin_mlat, vmag, dvmag, vdir, dvdir, chi2,
                            err_thres=100., mag_thres=100., titles=titles,
                            ylabel='Apex MLAT', yticklabels=binmlat, clim=clim, cmap=cmap,
                            filename=os.path.join(self.plotsavedir,vmag_fname))

            titles = ['E mag. (mV/m)', 'E mag err (mV/m)', 'E dir (deg)', 'E dir err (deg)', '']
            clim = [[0.,75.],[0., 15.],[-180., 180.],[0., 35.]]
            cmap = ['viridis', 'turbo', 'twilight', 'turbo']

            summary_plots.plot_magnitude(times, self.bin_mlat, emag, demag, edir, dedir, chi2,
                            err_thres=5., mag_thres=5., titles=titles,
                            ylabel='Apex MLAT', yticklabels=binmlat, clim=clim, cmap=cmap,
                            filename=os.path.join(self.plotsavedir,emag_fname))
