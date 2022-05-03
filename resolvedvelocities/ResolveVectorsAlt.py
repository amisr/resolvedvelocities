#!/usr/bin/env python
"""


Created: 4 July 2018

@author: Ashton S. Reimer


0.......10........20........30........40........50........60........70........80
"""

import os
import sys
import tables
import numpy as np
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import socket
import getpass
import platform
from argparse import ArgumentParser, RawDescriptionHelpFormatter


from .DataHandler import FittedVelocityDataHandler
from .plot import summary_plots
from .utils import *
import resolvedvelocities as rv

config_file_help = """Calculate 2D resolved plasma drift velocity and electric
field vectors from the LoS measurments in a fitted AMISR file in altitude bins.

Requires a configuration file containing the following example format:\n""" + get_example_config()



class ResolveVectorsAlt(object):
    def __init__(self,configfile):

        # read the config file
        self.configfile = configfile
        self.read_config(configfile)

        print(self.datafile)

        self.datahandler = FittedVelocityDataHandler(self.datafile)
        self.datahandler.load_data(self.use_beams)
        self.datahandler.filter(chi2=self.chi2lim, ne=self.nelim, mlat=self.mlatlim, fitcode=self.goodfitcode, chirp=self.chirp)

        if self.integration_time:
            self.integration_periods = get_integration_periods(self.datahandler.utime, self.integration_time)
        else:
            self.integration_periods = self.datahandler.utime

        self.altitude_bins = self.get_bins()



    def run(self):

        self.transform()
        self.compute_vector_velocity()
        self.compute_geodetic_velocity()
        self.save_output()
        if self.plotsavedir:
            self.create_plots()


    def read_config(self, config_file):

        # read config file
        config = configparser.ConfigParser(converters={'list':parse_list})
        config.read(config_file)

        # input/output file paths
        self.datafile = config.get('FILEIO', 'DATAFILE')
        self.output_path = config.get('FILEIO', 'OUTPUT_PATH')
        self.output_name = config.get('FILEIO', 'OUTPUT_NAME')

        # general options
        self.chirp = config.getfloat('OPTIONS', 'CHIRP') if config.has_option('OPTIONS', 'CHIRP') else 0.
        self.nelim = config.getlist('OPTIONS', 'NELIM') if config.has_option('OPTIONS', 'NELIM') else None
        self.chi2lim = config.getlist('OPTIONS', 'CHI2LIM') if config.has_options('OPTIONS', 'CHI2LIM') else None
        self.goodfitcode = config.getlist('OPTIONS', 'GOODFITCODE') if config.has_option('OPTIONS', 'GOODFITCODE') else None
        self.aprior_covar = config.getlist('OPTIONS', 'COVAR')
        self.integration_time = config.getfloat('OPTIONS', 'INTTIME') if config.has_option('OPTIONS', 'INTTIME') else None
        self.use_beams = config.getlist('OPTIONS', 'USE_BEAMS') if config.has_option('OPTIONS', 'USE_BEAMS') else None

        # altitude-resolved vector velocities specific options
        self.altitude_bins_def = config.get('VVELSALT', 'ALTBIN')
        self.mlatlim = config.getlist('VVELSALT', 'MLATLIM') if config.has_option('VVELS_ALT', 'MLATLIM') else None

        # plotting
        self.plotsavedir = config.get('PLOTTING', 'PLOTSAVEDIR') if config.has_option('PLOTTING', 'PLOTSAVEDIR') else None



    # bin forming function
    def get_bins(self):

        self.bins = np.empty((0,2))
        self.bin_centers = np.empty((0,))

        bins_groups = self.altitude_bins_def.split(';')
        for i,group in enumerate(bins_groups):
            start, stop, step, stride = [float(i) for i in group.split(',')]

            bin_centers0 = np.arange(start, stop, step)
            bins0 = np.array([bin_centers0-stride/2., bin_centers0+stride/2.]).T

            self.bins = np.append(self.bins, bins0, axis=0)
            self.bin_centers = np.append(self.bin_centers, bin_centers0)




    def rotation_matrices(self, lat, lon):
        lam = lat*np.pi/180.
        phi = lon*np.pi/180.

        R_enu2uvw = np.array([[-np.sin(phi), -np.sin(lam)*np.cos(phi), np.cos(lam)*np.cos(phi)],
                    [np.cos(phi), -np.sin(lam)*np.sin(phi), np.cos(lam)*np.sin(phi)],
                    [np.cos(lam), np.zeros(lam.shape), np.sin(lam)]])

        R_uvw2enu = np.array([[-np.sin(phi), np.cos(phi), np.zeros(lam.shape)],
                   [-np.sin(lam)*np.cos(phi), -np.sin(lam)*np.sin(phi), np.cos(lam)],
                   [np.cos(lam)*np.cos(phi), np.cos(lam)*np.sin(phi), np.sin(lam)]])

        if np.ndim(lat)>0:
            R_enu2uvw = R_enu2uvw.transpose((2,0,1))
            R_uvw2enu = R_uvw2enu.transpose((2,0,1))

        return R_enu2uvw, R_uvw2enu

    def transform(self):

        R_enu2uvw, R_uvw2enu = self.rotation_matrices(self.datahandler.lat, self.datahandler.lon)

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([self.datahandler.ke, self.datahandler.kn, self.datahandler.kz]).T
        self.A = np.einsum('mij,...mj->...mi', R_enu2uvw, kvec)



    # For each integration period, bin the los velocity data by altitude
    # and then do the vvels on it
    # Code doesn't write anything to disk, uses RAM.
    def compute_vector_velocity(self):

        num_integrations = len(self.integration_periods)
        num_bins = len(self.bins)

        self.Velocity = np.full((num_integrations,num_bins,3), np.nan)
        self.VelocityCovariance = np.full((num_integrations,num_bins,3,3), np.nan)
        self.ChiSquared = np.full((num_integrations,num_bins), np.nan)
        self.NumPoints = np.full((num_integrations,num_bins), np.nan)


        for i,integration_period in enumerate(self.integration_periods):

            tidx = self.datahandler.get_record_indices(integration_period[0],integration_period[1])

            # do vvels for each altitude bin
            for k, (bin_start, bin_end) in enumerate(self.bins):

                aidx = np.where((self.datahandler.alt >= bin_start) & (self.datahandler.alt <= bin_end))[0]
                if len(aidx) == 0:
                    continue
                vlos = self.datahandler.vlos[tidx,aidx[:,np.newaxis]].flatten()
                dvlos = self.datahandler.dvlos[tidx,aidx[:,np.newaxis]].flatten()
                A = np.repeat(self.A[aidx], len(tidx), axis=0)

                V, SigV, chi2, N = vvels(vlos, dvlos, A, self.aprior_covar, minnumpoints=3)

                # now transform from local uvw to enu
                self.Velocity[i,k,:] = V
                self.VelocityCovariance[i,k,:] = SigV
                self.ChiSquared[i,k] = chi2
                self.NumPoints[i,k] = N


    def compute_geodetic_velocity(self):

        R_enu2uvw, R_uvw2enu = self.rotation_matrices(self.datahandler.site[0], self.datahandler.site[1])


        # now transform from local uvw to enu
        self.Velocity = np.einsum('ij,...mj->...mi', R_uvw2enu, self.Velocity)
        self.VelocityCovariance = np.einsum('ij,...mjk,lk->...mil', R_uvw2enu, self.VelocityCovariance, R_uvw2enu)

        # calculate vector magnitude and direction
        north = np.tile(np.array([0,1,0]), (len(self.bin_centers),1))
        self.Vmag, self.Vmag_err, self.Vdir, self.Vdir_err = magnitude_direction(self.Velocity, self.VelocityCovariance, north)




    def save_output(self):

        # save output file
        os.makedirs(os.path.abspath(self.output_path),exist_ok=True)
        outfile = os.path.join(self.output_path, self.output_name)
        FILTERS = tables.Filters(complib='zlib', complevel=1)
        with tables.open_file(outfile, mode='w',filters=FILTERS) as outfile:

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


            outfile.create_group('/', 'VectorVels')
            save_carray(outfile, '/VectorVels/AltitudeBins', self.bins, {'TITLE':'Altitude Bins'})
            save_carray(outfile, '/VectorVels/Altitudes', self.bin_centers, {'TITLE':'Altitude of Center of Bins'})
            save_carray(outfile, '/VectorVels/Velocity', self.Velocity, {'TITLE':'Plasma Drift Velocity', 'Size':'Nrecords x Nalts x 3 (East, North, Up)', 'Units':'m/s'})
            save_carray(outfile, '/VectorVels/CovarianceV', self.VelocityCovariance, {'TITLE':'Velocity Covariance Matrix', 'Size':'Nrecords x Nalts x 3 x 3', 'Units':'(m/s)^2'})
            save_carray(outfile, '/VectorVels/Vmag', self.Vmag, {'TITLE':'Velocity Magnitude', 'Size':'Nrecords x Nalts', 'Units':'m/s'})
            save_carray(outfile, '/VectorVels/errVmag', self.Vmag_err, {'TITLE':'Velocity Magnitude Error', 'Size':'Nrecords x Nalts', 'Units':'m/s'})
            save_carray(outfile, '/VectorVels/Vdir', self.Vdir, {'TITLE':'Velocity Direction Angle East of North', 'Size':'Nrecords x Nalts', 'Units':'Degrees'})
            save_carray(outfile, '/VectorVels/errVdir', self.Vdir_err, {'TITLE':'Velocity Direction Error', 'Size':'Nrecords x Nalts', 'Units':'Degrees'})


            outfile.create_group('/', 'ProcessingParams')

            save_carray(outfile, '/ProcessingParams/Chi2', self.ChiSquared, {'TITLE':'Reduced Chi-Squared', 'Size':'Nrecords x Nbins'})
            save_carray(outfile, '/ProcessingParams/NumPoints', self.NumPoints, {'TITLE':'Number of input data points used to estimate the vector', 'Size':'Nrecords x Nbins'})

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


    def create_plots(self):

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
                vcom_fname = 'alt_velocity_components.png'
                vmag_fname = 'alt_velocity_magnitudes.png'
            else:
                vcom_fname = 'alt_velocity_components_{}.png'.format(t)
                vmag_fname = 'alt_velocity_magnitudes_{}.png'.format(t)

            # make vector plots
            times = self.integration_periods[start_ind:end_ind,:]
            binmlat = ['{:.2f}'.format(ml) for ml in self.bin_centers]

            vels = self.Velocity[start_ind:end_ind,:]
            covvels = self.VelocityCovariance[start_ind:end_ind,:]

            summary_plots.plot_components(times, self.bin_centers, vels, covvels,
                            titles=[r'$V_E$ (m/s)',r'$V_N$ (m/s)',r'$V_U$ (m/s)'],
                            ylabel='Altitude (km)', yticklabels=binmlat,
                            clim=[[-1500.,1500.], [0.,350.]], cmap=['coolwarm', 'turbo'],
                            filename=os.path.join(self.plotsavedir,vcom_fname), scale_factors=[1,1,1])



            # make magnitude plots
            vmag = self.Vmag[start_ind:end_ind,:]
            dvmag = self.Vmag_err[start_ind:end_ind,:]
            vdir = self.Vdir[start_ind:end_ind,:]
            dvdir = self.Vdir_err[start_ind:end_ind,:]
            chi2 = self.ChiSquared[start_ind:end_ind,:]

            titles = ['V mag. (m/s)', 'V mag. err. (m/s)', 'V dir. (deg)', 'V dir. err. (deg)', '']
            clim = [[0.,1500.],[0., 350.],[-180., 180.],[0., 35.]]
            cmap = ['viridis', 'turbo', 'twilight', 'turbo']

            summary_plots.plot_magnitude(times, self.bin_centers, vmag, dvmag, vdir, dvdir, chi2,
                            err_thres=100., mag_thres=100., titles=titles,
                            ylabel='Altitude (km)', yticklabels=binmlat, clim=clim, cmap=cmap,
                            filename=os.path.join(self.plotsavedir,vmag_fname))


def main():
    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('config_file',help='A configuration file.')

    args = vars(parser.parse_args())
    vvelsalt = ResolveVectorsAlt(args['config_file'])
    vvelsalt.run()


if __name__ == '__main__':
    main()
