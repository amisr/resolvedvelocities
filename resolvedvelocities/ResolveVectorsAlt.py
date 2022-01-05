#!/usr/bin/env python
"""


Created: 4 July 2018

@author: Ashton S. Reimer


0.......10........20........30........40........50........60........70........80
"""
from __future__ import print_function
from __future__ import division

import re
import os
import sys
import glob
import tables
import shutil
import numpy as np
from datetime import datetime
import tempfile
try:
    import ConfigParser as configparser
except ImportError:
    import configparser
import socket
import getpass
import platform

import resolvedalts as ra
from resolvedvelocities.DataHandler import FittedVelocityDataHandler
import resolvedvelocities as rv

# output file definition
h5paths = [['/ProcessingParams','Processing Parameters'],
           ['/VectorVels','Velocity Vector Reconstructed From LOS'],
           ['/Site','Site Parameters'],
           ['/Time','Time Information'],
          ]

h5attribs = {'/VectorVels/AltitudeBins' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/Altitudes' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/Vest' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/covVest' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/errVest' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/Nmeas' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/VectorVels/chi2' : [('TITLE','Altitude'),('Unit','Meters'),('Description','Altitude assuming local flat Earth.')],
             '/ProcessingParams/ProcessingTimeStamp' : [('TITLE','Processing Time Stamp')],
             '/Site/Altitude' : [('TITLE','Altitude'),('Description','Altitude of site'),('Unit','Meters')],
             '/Site/Code' : [('TITLE','Site Code')],
             '/Site/Latitude' : [('TITLE','Latitude'),('Description','Latitude of site'),('Unit','Degrees North')],
             '/Site/Longitude' : [('TITLE','Longitude'),('Description','Longitude of site'),('Unit','Degrees East')],
             '/Site/Name' : [('TITLE','Name'),('Description','Site Name')],
             '/Time/Day' : [('TITLE','Day of Month'),('Size','Nrecords x 2 (Start and end of integration')],
             '/Time/Month' : [('TITLE','Month'),('Size','Nrecords x 2 (Start and end of integration')],
             '/Time/Year' : [('TITLE','Year'),('Size','Nrecords x 2 (Start and end of integration')],
             '/Time/doy' : [('TITLE','Day of Year'),('Size','Nrecords x 2 (Start and end of integration')],
             '/Time/UnixTime' : [('TITLE','Unix Time'),('Size','Nrecords x 2 (Start and end of integration'),('Unit','Seconds')],
            }


class ResolveVectorsAlt(object):
    def __init__(self,configfile):

        # read the config file
        print("Reading configuration file...")
        # self.configfile = configfile
        # self.config = self.parse_config()
        self.configfile = configfile
        self.read_config(configfile)

        # move to config file
        self.aprior_covar = np.array([3000.0**2,3000.0**2,100.0**2])

        # # find files using config information
        # print("Checking input files...")
        # self.filelists = self.get_files()
        #
        # # initialize the data handlers
        # print("Checking available data...")
        # self.num_pulsetypes = len(self.filelists)
        # # self.datahandlers = [FittedVelocityDataHandler(self.filelists[i]) for i in range(self.num_pulsetypes)]

        print(self.datafile)
        self.datahandler = FittedVelocityDataHandler(self.datafile)
        self.datahandler.read_data(self.use_beams)
        self.datahandler.filter(chi2=self.chi2lim, ne=self.nelim, alt=self.altlim, fitcode=self.goodfitcode, chirp=self.chirp)


        # set up the bins
        print("Forming altitude bins...")
        self.altitude_bins = self.get_bins()

        # find all unique times from the data handlers and then
        # determine the integration periods to calculate Ne_NoTr for
        print("Calculating integration periods...")
        # self.times = self.get_unique_times()
        self.integration_periods = self.get_integration_periods()

        # make sure output directory is available and if not create it
        print("Validating output directory...")
        # output_dir = self.config['output']['output_path']
        # if not os.path.exists(output_dir):
        if not os.path.exists(self.output_path):
            print("    Output directory doesn't exist!")
            print("    Attempting to create one...")
            os.makedirs(self.output_path)


    def read_config(self, config_file):

        # read config file
        config = configparser.ConfigParser(converters={'list':self.parse_list})
        config.read(config_file)

        # Possibly could done better with converters?  This may by python3 specific though.
        self.datafile = config.get('FILEIO', 'files')
        self.output_path = config.get('FILEIO', 'OUTPUT_PATH')
        self.output_name = config.get('FILEIO', 'OUTPUT_NAME')

        self.altitude_bins_def = config.get('VVELS_OPTIONS', 'altitude_bins')
        self.chirp = config.getfloat('VVELS_OPTIONS', 'CHIRP')
        self.aprior_covar = config.getlist('VVELS_OPTIONS', 'COVAR')
        self.altlim = config.getlist('VVELS_OPTIONS', 'ALTLIM')
        self.nelim = config.getlist('VVELS_OPTIONS', 'NELIM')
        self.chi2lim = config.getlist('VVELS_OPTIONS', 'CHI2LIM')
        self.goodfitcode = config.getlist('VVELS_OPTIONS', 'GOODFITCODE')
        self.integration_time = config.getfloat('VVELS_OPTIONS', 'INTTIME') if config.has_option('VVELS_OPTIONS', 'INTTIME') else None
        self.use_beams = config.getlist('VVELS_OPTIONS', 'USE_BEAMS') if config.has_option('VVELS_OPTIONS', 'USE_BEAMS') else None

        # self.outfilename = config.get('FILEIO', 'OUTFILENAME')
        # self.chirp = config.getfloat('CONFIG', 'CHIRP')
        # self.covar = [float(i) for i in config.get('CONFIG', 'COVAR').split(',')]
        # self.altlim = [float(i) for i in config.get('CONFIG', 'ALTLIM').split(',')]
        # self.nelim = [float(i) for i in config.get('CONFIG', 'NELIM').split(',')]
        # self.chi2lim = [float(i) for i in config.get('CONFIG', 'CHI2LIM').split(',')]
        # self.goodfitcode = [float(i) for i in config.get('CONFIG', 'GOODFITCODE').split(',')]
        # self.binvert = eval(config.get('CONFIG', 'BINVERT'))
        # # can probably change this parameter to a list of start, end, step similar to vvelsAlt
        # self.outalt = [float(i) for i in config.get('CONFIG', 'OUTALT').split(',')]
        # self.marprot = [float(i) for i in config.get('CONFIG', 'MARPROT').split(',')]
        #
        # # optional parameters
        # self.plotsavedir = config.get('PLOTTING', 'PLOTSAVEDIR') if config.has_option('PLOTTING', 'PLOTSAVEDIR') else None
        # self.upB_beamcode = config.getint('CONFIG', 'UPB_BEAMCODE') if config.has_option('CONFIG', 'UPB_BEAMCODE') else None
        # self.ionup = config.get('CONFIG', 'IONUP') if config.has_option('CONFIG', 'IONUP') else None
        # self.use_beams = [int(i) for i in config.get('CONFIG', 'USE_BEAMS').split(',')] if config.has_option('CONFIG', 'USE_BEAMS') else None
        # self.integration_time = config.getfloat('CONFIG', 'INTTIME') if config.has_option('CONFIG', 'INTTIME') else None
        # self.outfilepath = config.get('FILEIO', 'OUTFILEPATH') if config.has_option('FILEIO', 'OUTFILEPATH') else '.'

    def parse_list(self, s):
        return [float(i) for i in s.split(',')]

    # def parse_config(self):
    #     # Use ALLCAPS convection for config file
    #     required_sections_options = {'input': {'files': str},
    #                                  'output': {'output_name':str,
    #                                             'output_path': str,
    #                                            },
    #                                  'vvels_options': {'recs2integrate': int,
    #                                                    'altitude_bins':str},
    #                                 }
    #
    #     optional_sections_options = {'input': {'amb_path': str,
    #                                            'ksys_file': str
    #                                           },
    #                                 }
    #
    #
    #     # read the config file and convert to dictionary
    #     parser = ConfigParser.ConfigParser()
    #     parser.read(self.configfile)
    #     parsed_config = self.__config_to_dict_helper(parser)
    #
    #     # check the config file to make sure we have all required information
    #     # This checking + error statements handled by configparser?
    #     for section in required_sections_options.keys():
    #         if parsed_config.get(section,None) is None:
    #             msg = 'Required section: "%s" is missing from config.' % section
    #             raise AttributeError(msg)
    #         for option in required_sections_options[section].keys():
    #             if parsed_config[section].get(option,None) is None:
    #                 msg = 'Required option: "%s" is missing' % option
    #                 msg += ' from the "%s" section in the config.' % section
    #                 raise AttributeError(msg)
    #
    #             # potentially useful way to solve conversion/typing issues
    #             # convert the input config data to the required format
    #             type_func = required_sections_options[section][option]
    #             converted = type_func(parsed_config[section][option])
    #             parsed_config[section][option] = converted
    #
    #     # make sure optional options are formatted as required
    #     for section in optional_sections_options.keys():
    #         for option in optional_sections_options[section].keys():
    #             # convert the input config data to the required format
    #             type_func = optional_sections_options[section][option]
    #             try:
    #                 converted = type_func(parsed_config[section][option])
    #                 parsed_config[section][option] = converted
    #             except KeyError:
    #                 pass
    #
    #     return parsed_config


    def get_files(self):
        # we need to find all files that match the search strings
        # and check every input file path for them
        # This combines AC/LP files?
        # by_pulsetype = self.config['input']['files'].split(',')
        by_pulsetype = self.datafile.split(',')
        num_pulsetypes = len(by_pulsetype)

        filelists = [[] for x in range(num_pulsetypes)]
        for i in range(num_pulsetypes):
            files_for_pulsetype = by_pulsetype[i].split(':')
            # check that all files are accessible
            for f in files_for_pulsetype:
                if not os.path.exists(f):
                    raise FileNotFoundError('File does not exist: %s' % f)
            filelists[i].extend(sorted(files_for_pulsetype))

        num_pulsetypes = len(filelists)

        # calculate the total number of files
        num_files = 0
        for i in range(len(filelists)):
            num_files += len(filelists[i])

        filestr = "files" if num_files > 1 else "file"
        pulsetypestr = "pulsetypes" if num_pulsetypes > 1 else "pulsetype"
        print("Found %s %s for %s %s." % (num_files,filestr,num_pulsetypes,pulsetypestr))

        return filelists


    # bin forming function
    def get_bins(self):
        # bins_groups = self.config['vvels_options']['altitude_bins'].split(':')
        bins_groups = self.altitude_bins_def.split(':')
        for i,group in enumerate(bins_groups):
            opts = group.split(',')
            start  = float(opts[0])
            stop   = float(opts[1])
            step   = float(opts[2])
            stride = float(opts[3])
            num_alts = np.floor((stop-start)/step)+1
            temp_bin_centers = np.arange(start,stop,step=step)
            temp_bins = np.repeat(temp_bin_centers[:,np.newaxis],2,axis=1)
            temp_bins[:,0] -= stride / 2.
            temp_bins[:,1] += stride / 2.
            if i == 0:
                bins = temp_bins
                bin_centers = temp_bin_centers
            else:
                bins = np.concatenate((bins,temp_bins),axis=0)
                bin_centers = np.concatenate((bin_centers,temp_bin_centers),axis=0)

        self.bins = bins
        self.bin_centers = bin_centers


    def get_unique_times(self):
        all_times = list()
        for i in range(self.num_pulsetypes):
            all_times.extend(list(self.datahandlers[i].times))
        all_times = [tuple(x) for x in all_times]
        unique_times = np.array(sorted(list(set(list(all_times)))))

        # now detect time pairs that have 0 difference in start or end time
        # sometimes fitted files don't have exactly the same time windows...
        if self.num_pulsetypes > 1:
            diffs = np.diff(unique_times,axis=0)   # diff the start and end times
            diffs = np.array([[x[0].total_seconds(),x[1].total_seconds()] for x in diffs])
            inds = np.where(~((diffs[:,0] == 0) | (diffs[:,1] == 0)))[0]  # exclude times where start or end diffs are 0
            unique_times = unique_times[inds,:]

        return unique_times


    def get_integration_periods(self):

        if not self.integration_time:
            return self.datahandler.data['time']

        integration_periods = list()
        start_time = None
        # integration_time = self.config['vvels_options']['recs2integrate']
        integration_time = self.integration_time
        num_times = len(self.datahandler.data['time'])
        for i,time_pair in enumerate(self.datahandler.data['time']):
        # integration_time = self.integration_time
        # num_times = len(self.times)
        # for i,time_pair in enumerate(self.times):
            temp_start_time, temp_end_time = time_pair
            if start_time is None:
                start_time = temp_start_time
            time_diff = (temp_end_time - start_time).total_seconds()

            if time_diff >= integration_time:
                integration_periods.append([start_time,temp_end_time])
                start_time = None
                continue

            # Add an integration period for when we are at the end of the files
            # but we haven't reached the requested integration time
            if (i == num_times -1):
                integration_periods.append([start_time,temp_end_time])

        return np.array(integration_periods)


    @staticmethod
    # configparser objects can basically be treated as a dictionary
    # this may be outdated?
    def __config_to_dict_helper(configparserclass):
        # based on https://gist.github.com/amitsaha/3065184
        # converts a config parser object to a dictionary
        config = dict()
        defaults = configparserclass.defaults()
        sections = configparserclass.sections()

        temp = dict()
        for key in defaults:
            temp[key] = defaults[key]
        config['default'] = temp
        default_options = temp.keys()

        for section in sections:
            opts = configparserclass.options(section)
            options = [x for x in opts if not x in default_options]
            temp = dict()
            for option in options:
                temp[option] = configparserclass.get(section,option)
            config[section.lower()] = temp

        return config

    def transform(self):
        # lon = data['site']['lon']*np.pi/180.
        # neu_to_xyz = np.matrix([[np.cos(lon), -np.sin(lon), 0],[np.sin(lon),np.cos(lon),0],[0,0,1]])
        # xyz_to_neu = np.linalg.inv(neu_to_xyz)
        #
        # kvec_xyz = data['kvec'] #np.array((neu_to_xyz * np.matrix(data['kvec']).T).T)

        lam = self.datahandler.data['lat']*np.pi/180.
        phi = self.datahandler.data['lon']*np.pi/180.

        R_enu2uvw = np.array([[-np.sin(phi), -np.sin(lam)*np.cos(phi), np.cos(lam)*np.cos(phi)],
                    [np.cos(phi), -np.sin(lam)*np.sin(phi), np.cos(lam)*np.sin(phi)],
                    [np.cos(lam), np.zeros(lam.shape), np.sin(lam)]]).transpose((2,0,1))

        R_uvw2enu = np.array([[-np.sin(phi), np.cos(phi), np.zeros(lam.shape)],
                   [-np.sin(lam)*np.cos(phi), -np.sin(lam)*np.sin(phi), np.cos(lam)],
                   [np.cos(lam)*np.cos(phi), np.cos(lam)*np.sin(phi), np.sin(lam)]]).transpose((2,0,1))


        # glat = self.datahandler.data['lat']
        # glon = self.datahandler.data['lon']
        # galt = self.datahandler.data['alt'] / 1000.

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([self.datahandler.data['ke'], self.datahandler.data['kn'], self.datahandler.data['kz']]).T
        self.A = np.einsum('mij,...mj->...mi', R_enu2uvw, kvec)

        # u, v, w = enu2uvw(self.datahandler.data['ke'], self.datahandler.data['kn'], self.datahandler.data['kz'], glat, glon)
        #
        # self.A = np.array([u, v, w]).T


    # For each integration period, bin the los velocity data by altitude
    # and then do the vvels on it
    # Code doesn't write anything to disk, uses RAM.
    def do_vvels(self):

        num_integrations = len(self.integration_periods)
        num_bins = len(self.bins)

        vvels = np.zeros((num_integrations,num_bins,3))*np.nan
        chi2s = np.zeros((num_integrations,num_bins))*np.nan
        num_points = np.zeros((num_integrations,num_bins))*np.nan
        vvels_cov = np.zeros((num_integrations,num_bins,3,3))*np.nan
        for i,integration_period in enumerate(self.integration_periods):
            print("Integration period %s/%s" % (str(i+1),str(num_integrations)))

            kvecs = list()
            vels = list()
            evels = list()

            # for j,datahandler in enumerate(self.datahandlers):
            data, _ = self.datahandler.get_records(integration_period[0],integration_period[1])

                # can probably use pymap3d here
                # to prevent distortion due to lat/lon grid, we want to convert
                # the local North, East, Up kvectors to local x, y, z. Basically
                # just a rotation to a polar grid, so build the transformation matrix
                # if i == 0 and j == 0:
                #     lon = data['site']['lon']*np.pi/180.
                #     neu_to_xyz = np.matrix([[np.cos(lon), -np.sin(lon), 0],[np.sin(lon),np.cos(lon),0],[0,0,1]])
                #     xyz_to_neu = np.linalg.inv(neu_to_xyz)
                #
                # # # move filtering to datahandler?
                # # # definitely move hardcoded limits to config file
                # # # filter data by Ne < 10**9 and by chi2 < 10
                # # bad_ne_inds = np.where(data['ne'] < 1.0e9)
                # # bad_chi2_inds = np.where(data['chi2'] > 10)
                # # bad_mlat_inds = np.where(data['mlat'] > 69.0)
                # # data['vel'][bad_ne_inds] = np.nan
                # # data['evel'][bad_ne_inds] = np.nan
                # # data['vel'][bad_chi2_inds] = np.nan
                # # data['evel'][bad_chi2_inds] = np.nan
                # # data['vel'][:,bad_mlat_inds] = np.nan
                # # data['evel'][:,bad_mlat_inds] = np.nan
                #
                # kvec_xyz = data['kvec'] #np.array((neu_to_xyz * np.matrix(data['kvec']).T).T)

                # # Lot of array manipulation here?
                # # now reshape into columns
                # num_times, num_beams, num_alts = data['vel'].shape
                #
                # # flatten arrays into a column of measurements
                # temp_vel = np.swapaxes(data['vel'],0,1).flatten()
                # temp_evel = np.swapaxes(data['evel'],0,1).flatten()
                # temp_alt = np.repeat(data['altitude'][:,np.newaxis,:],num_times,axis=1).flatten()
                # # kvector can't be flattened, we need an array of
                # # num_measurements x 3
                # temp_kvec = list()
                # for k in range(num_beams):
                #     temp_k = [kvec_xyz[k,:].tolist() for x in range(num_times*num_alts)]
                #     temp_kvec.extend(temp_k)
                # temp_kvec = np.array(temp_kvec)
                #
                # if j == 0:
                #     vels  = temp_vel
                #     evels = temp_evel
                #     alts  = temp_alt
                #     kvecs = temp_kvec
                # else:
                #     vels  = np.concatenate((vels,temp_vel),axis=0)
                #     evels = np.concatenate((evels,temp_evel),axis=0)
                #     alts  = np.concatenate((alts,temp_alt),axis=0)
                #     kvecs = np.concatenate((kvecs,temp_kvec),axis=0)

            # do vvels for each altitude bin
            for k, (bin_start, bin_end) in enumerate(self.bins):

                alt_inds = np.where((self.datahandler.data['alt'] >= bin_start) & (self.datahandler.data['alt'] <= bin_end))[0]
                if len(alt_inds) == 0:
                    continue
                input_vlos = data['vels'][alt_inds].T
                input_dvlos = data['evels'][alt_inds]
                # A = kvecs[alt_inds,:]

                v, sigv, chi2, N = self.vvels(input_vlos, input_dvlos, self.A, self.aprior_covar, minnumpoints=3)

                # now transform from local x,y,z to N,E,U
                # np.array((neu_to_xyz * np.matrix(data['kvec']).T).T)
                # THIS CONVERSION OF COVARIANCE MATRIX (ECEF->GD) YOU HAD TO DO FOR AGU POSTER
                # LOOK UP FROM GRADIENT CALCULATION WORK
                vvels[i,k,:] = np.einsum('mij,...mj->...mi', R_uvw2enu, v)
                vvels_cov[i,k,:] = np.einsum('mij,...mjk,mlk->...mil',R_uvw2enu, sigv, R_uvw2enu)

                # vvels[i,k,:]     = (xyz_to_neu * np.matrix(v).T).T
                # vvels_cov[i,k,:] = xyz_to_neu * sigv * xyz_to_neu.T
                chi2s[i,k]       = chi2
                num_points[i,k]  = N

        data = dict()
        data['vvels'] = vvels
        data['chi2s'] = chi2s
        data['num_points'] = num_points
        data['vvels_cov'] = vvels_cov

        return data


# same vvels function?  Move to a utils file?
    def vvels(self, vlos, dvlos, A, cov, minnumpoints=1):
        # implimentation of Heinselman and Nicolls 2008 vector velocity algorithm

        # remove nan data points
        finite_inds = np.where(np.isfinite(vlos))[0]

        # if there are too few points for a valid reconstruction, set output to NAN
        if finite_inds.size < minnumpoints:
            V = np.full(3,np.nan)
            sigV = np.full((3,3),np.nan)
            return V, sigV, np.nan, 0.

        vlos = vlos[finite_inds]
        dvlos = dvlos[finite_inds]
        A = np.matrix(A[finite_inds,:])

        sigmaE = np.matrix(np.diag(dvlos**2))
        sigmaV = np.matrix(np.diag(cov))

        try:
            X = sigmaV*A.T*np.linalg.inv(A*sigmaV*A.T + sigmaE)
            V = X*np.matrix(vlos).T     # calculate velocity estimate (Heinselman 2008 eqn 12)
            sigV = sigmaV - X*A*sigmaV  # calculate covariance of velocity estimate (Heinselman 2008 eqn 13)
            temp = np.array(A*V - np.matrix(vlos).T)
            chi2 = np.sum(temp**2/dvlos**2) / finite_inds.size
            V = np.array(V.T)
        except np.linalg.LinAlgError:
            V = np.full(3,np.nan)
            sigV = np.full((3,3),np.nan)
            chi2 = np.nan

        return V, sigV, chi2, finite_inds.size


    def get_site(self):
        site = dict()
        fname = self.datahandler.filelist
        with tables.open_file(fname,'r') as h5:
            site['altitude']  = h5.root.Site.Altitude.read()
            site['code']      = h5.root.Site.Code.read()
            site['latitude']  = h5.root.Site.Latitude.read()
            site['longitude'] = h5.root.Site.Longitude.read()
            site['name']      = h5.root.Site.Name.read()

        return site


    # def get_time(self):
    #     epoch = datetime(1970,1,1)
    #     time_shape = self.integration_periods.shape
    #     time = dict()
    #     keys = ['day','month','year','doy','unixtime']
    #     for key in keys:
    #         time[key] = np.zeros(time_shape,dtype=int)
    #
    #     for i,pair in enumerate(self.integration_periods):
    #         time['day'][i,:]   = np.array([pair[0].day,pair[1].day])
    #         time['month'][i,:] = np.array([pair[0].month,pair[1].month])
    #         time['year'][i,:]  = np.array([pair[0].year,pair[1].year])
    #         time['doy'][i,:]   = np.array([pair[0].timetuple().tm_yday,pair[1].timetuple().tm_yday])
    #         diff_pair = [(pair[0]-epoch).total_seconds(),(pair[1]-epoch).total_seconds()]
    #         time['unixtime'][i,:] = diff_pair
    #
    #     return time

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


    def save_output(self, output, errvvels):
        # TODO: come up with a better way to manage all this

        # save output file
        outfile = os.path.join(self.output_path, self.output_name)
        FILTERS = tables.Filters(complib='zlib', complevel=1)
        with tables.open_file(outfile, mode='w',filters=FILTERS) as outfile:

            # copy some groups directly from fitted input file
            with tables.open_file(self.datafile, mode='r') as infile:
                outfile.copy_children(infile.get_node('/Site'), outfile.create_group('/','Site'))
                if not self.integration_time:
                    outfile.copy_children(infile.get_node('/Time'), outfile.create_group('/','Time'))
                else:
                    outfile.create_group('/','Time')
                    year, month, day, doy, dtime, mlt = self.create_time_arrays()

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
            save_carray(outfile, '/VectorVels/Vest', output['vvels'], {'TITLE':'Velocity'})
            save_carray(outfile, '/VectorVels/covVest', output['vvels_cov'], {'TITLE':'Velocity Covariance Matrix'})
            save_carray(outfile, '/VectorVels/errVest', errvvels, {'TITLE':'Velocity Errors'})
            save_carray(outfile, '/VectorVels/Nmeas', output['num_points'], {'TITLE':'Number of Measurements'})
            save_carray(outfile, '/VectorVels/chi', output['chi2s'], {'TITLE':'Chi Squared'})


            # outfile.create_group('/', 'Magnetic')
            #
            # save_carray(outfile, '/Magnetic/Latitue', self.bin_alat, {'TITLE':'Magnetic Latitude', 'Size':'Nbins'})
            # save_carray(outfile, '/Magnetic/Longitude', self.bin_alon, {'TITLE':'Magnetic Longitude', 'Size':'Nbins'})
            # save_carray(outfile, '/Magnetic/Velocity', self.Velocity, {'TITLE':'Plasma Drift Velocity', 'Size':'Nrecords x Nbins x 3 (Ve1, Ve2, Ve3)', 'Units':'m/s'})
            # save_carray(outfile, '/Magnetic/CovarianceV', self.VelocityCovariance, {'TITLE':'Velocity Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(m/s)^2'})
            # save_carray(outfile, '/Magnetic/ElectricField', self.ElectricField, {'TITLE':'Convection Electric Field', 'Size':'Nrecords x Nbins x 3 (Ed1, Ed2, Ed3)', 'Units':'V/m'})
            # save_carray(outfile, '/Magnetic/CovarianceV', self.ElectricFieldCovariance, {'TITLE':'Electric Field Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(V/m)^2'})
            #
            #
            #
            # outfile.create_group('/', 'Geodetic')
            #
            # save_carray(outfile, '/Geodetic/Latitue', self.bin_alat, {'TITLE':'Geographic Latitude', 'Size':'Nalts x Nbins'})
            # save_carray(outfile, '/Geodetic/Longitude', self.bin_alon, {'TITLE':'Geographic Longitude', 'Size':'Nalts x Nbins'})
            # save_carray(outfile, '/Geodetic/Altitde', self.bin_alon, {'TITLE':'Geographic Altitude', 'Size':'Nalts x Nbins', 'Units':'km'})
            # save_carray(outfile, '/Geodetic/Velocity', self.Velocity_gd, {'TITLE':'Plasma Drift Velocity', 'Size':'Nrecords x Nalts x Nbins x 3 (East, North, Up)', 'Units':'m/s'})
            # save_carray(outfile, '/Geodetic/CovarianceV', self.VelocityCovariance_gd, {'TITLE':'Velocity Covariance Matrix', 'Size':'Nrecords x Nalts x Nbins x 3 x 3', 'Units':'(m/s)^2'})
            # save_carray(outfile, '/Geodetic/Vmag', self.Vgd_mag, {'TITLE':'Velocity Magnitude', 'Size':'Nrecords x Nalts x Nbins', 'Units':'m/s'})
            # save_carray(outfile, '/Geodetic/errVmag', self.Vgd_mag_err, {'TITLE':'Velocity Magnitude Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'m/s'})
            # save_carray(outfile, '/Geodetic/Vdir', self.Vgd_dir, {'TITLE':'Velocity Direction Angle East of North Magnetic Meridian (-e2)', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            # save_carray(outfile, '/Geodetic/errVdir', self.Vgd_dir_err, {'TITLE':'Velocity Direction Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            # save_carray(outfile, '/Geodetic/ElectricField', self.ElectricField_gd, {'TITLE':'Convection Electric Field', 'Size':'Nrecords x Nbins x 3 (East, North, Up)', 'Units':'V/m'})
            # save_carray(outfile, '/Geodetic/CovarianceE', self.ElectricFieldCovariance_gd, {'TITLE':'Electric Field Covariance Matrix', 'Size':'Nrecords x Nbins x 3 x 3', 'Units':'(V/m)^2'})
            # save_carray(outfile, '/Geodetic/Emag', self.Egd_mag, {'TITLE':'Electric Field Magnitude', 'Size':'Nrecords x Nalts x Nbins', 'Units':'V/m'})
            # save_carray(outfile, '/Geodetic/errEmag', self.Egd_mag_err, {'TITLE':'Electric Field Magnitude Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'V/m'})
            # save_carray(outfile, '/Geodetic/Edir', self.Egd_dir, {'TITLE':'Electric Field Direction Angle East of North Magnetic Meridian (-e2)', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})
            # save_carray(outfile, '/Geodetic/errEdir', self.Egd_dir_err, {'TITLE':'Electric Field Direction Error', 'Size':'Nrecords x Nalts x Nbins', 'Units':'Degrees'})


            # site = self.get_site()
            # outfile.create_group('/', 'Site')
            #
            # outfile.create_array('/Site', 'Altitude', site['altitude'])
            # h5.create_array('/Site', 'Code', site['code'])
            # h5.create_array('/Site', 'Latitude', site['latitude'])
            # h5.create_array('/Site', 'Longitude', site['longitude'])
            # h5.create_array('/Site', 'Name', site['name'])

            outfile.create_group('/', 'ProcessingParams')

            # RE-ADD THESE IN THE FUTURE
            # save_carray(outfile, '/ProcessingParams/Chi2', self.ChiSquared, {'TITLE':'Reduced Chi-Squared', 'Size':'Nrecords x Nbins'})
            # save_carray(outfile, '/ProcessingParams/NumPoints', self.NumPoints, {'TITLE':'Number of input data points used to estimate the vector', 'Size':'Nrecords x Nbins'})

            # outfile.create_group('/ProcessingParams', 'Apexpy')
            #
            # outfile.create_array('/ProcessingParams/Apexpy', 'Year', self.marp.year)
            # outfile.set_node_attr('/ProcessingParams/Apexpy/Year', 'TITLE', 'Decimal Year used for IGRF Model')
            #
            # outfile.create_array('/ProcessingParams/Apexpy','RefHeight',self.marp.refh)
            # outfile.set_node_attr('/ProcessingParams/Apexpy/RefHeight', 'TITLE', 'Reference height used for Apex coordinates')
            # outfile.set_node_attr('/ProcessingParams/Apexpy/RefHeight', 'Units', 'km')
            #
            # outfile.create_array('/ProcessingParams/Apexpy','Version',apexpy.__version__.encode('utf-8'))
            # outfile.set_node_attr('/ProcessingParams/Apexpy/Version', 'TITLE', 'Apexpy version used.')

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



# move to run_resolve_vectors.py?
    def run(self):

        # add overwrite flag to config file?
        # First check if output file is able to be created
        temp_file = tempfile.mktemp()
        # output_file = os.path.join(self.config['output']['output_path'],self.config['output']['output_name'])
        output_file = os.path.join(self.output_path,self.output_name)

        # Run the calculator and write output to a file
        output = self.do_vvels()
        i = np.array([0,1,2])
        errvvels = np.sqrt(np.squeeze(output['vvels_cov'][:,:,i,i]))

        # # get Site information
        # site = self.get_site()
        # # get Time information
        # time = self.get_time()

        self.save_output(output, errvvels)

        # # Get current date and time
        # date = datetime.utcnow()
        # processing_time = date.strftime("%a, %d %b %Y %H:%M:%S +0000")
        #
        # # move this to an independent function
        # # add more attribute information to the hdf5 output
        # # Write the output
        # # set up the output file
        # print("Writing data to file...")
        # with tables.open_file(temp_file,'w') as h5:
        #     for h5path in h5paths:
        #         group_path, group_name = os.path.split(h5path[0])
        #         h5.create_group(group_path,group_name,title=h5path[1],createparents=True)
        #
        #     node_path = '/ProcessingParams'
        #     h5.create_array(node_path,'ProcessingTimeStamp',np.array(processing_time),createparents=True)
        #
        #     node_path = '/Site'
        #     h5.create_array(node_path,'Altitude',site['altitude'],createparents=True)
        #     h5.create_array(node_path,'Code',site['code'],createparents=True)
        #     h5.create_array(node_path,'Latitude',site['latitude'],createparents=True)
        #     h5.create_array(node_path,'Longitude',site['longitude'],createparents=True)
        #     h5.create_array(node_path,'Name',site['name'],createparents=True)
        #
        #     node_path = '/VectorVels'
        #     h5.create_array(node_path,'AltitudeBins',self.bins,createparents=True)
        #     h5.create_array(node_path,'Altitudes',self.bin_centers,createparents=True)
        #     h5.create_array(node_path,'Vest',output['vvels'],createparents=True)
        #     h5.create_array(node_path,'covVest',output['vvels_cov'],createparents=True)
        #     h5.create_array(node_path,'errVest',errvvels,createparents=True)
        #     h5.create_array(node_path,'Nmeas',output['num_points'],createparents=True)
        #     h5.create_array(node_path,'chi2',output['chi2s'],createparents=True)
        #
        #     node_path = '/Time'
        #     h5.create_array(node_path,'Day',time['day'],createparents=True)
        #     h5.create_array(node_path,'Month',time['month'],createparents=True)
        #     h5.create_array(node_path,'Year',time['year'],createparents=True)
        #     h5.create_array(node_path,'doy',time['doy'],createparents=True)
        #     h5.create_array(node_path,'UnixTime',time['unixtime'],createparents=True)
        #
        #
        # # Add configuration information
        # print("Adding configuration information...")
        # inputfiles = [x.filelist for x in self.datahandlers]
        # self.write_config_info(temp_file,inputfiles)
        #
        # with tables.open_file(temp_file,'r+') as h5:
        #     for key in h5attribs.keys():
        #         for attr in h5attribs[key]:
        #             h5.set_node_attr(key,attr[0],attr[1])
        #
        # # Add repacking to vvelsLat
        # # repack the file with compression
        # print("Repacking the file with compression...")
        # repack = repackh5(temp_file,output_file)
        # repack.repack()
        # # remove the temporary file
        # print("Cleaning up...")
        # os.remove(temp_file)

        print("Plotting...")
        plotter = ra.ResolvedAltPlotter(output_file)
        plotter.make_plot()

        print("Done!")



    def write_config_info(self,h5name,raw_files):
        import platform
        import getpass

        # Computer information:
        PythonVersion   = platform.python_version()
        Type            = platform.machine()
        System          = "%s %s %s" % (platform.system(),platform.release(),platform.version())
        User            = getpass.getuser()
        Hostname        = platform.node()
        if len(Hostname) == 0:
            import socket
            Hostname = socket.gethostname()

        # Fitter configuration information
        Version = ra.__version__  # Eventually replaced by self.__version__

        # Get the config file used
        cf = self.configfile
        Path = os.path.dirname(os.path.abspath(cf))
        Name = os.path.basename(cf)

        with open(cf,'r') as f:
            Contents = "".join(f.readlines())

        # Record the raw files used
        # Make a string listing all the files
        RawFiles = ''
        for i,files in enumerate(raw_files):
            temp = "\n".join(files)
            if i != 0:
                RawFiles += '\n'
            RawFiles += temp

        # Record the directory where fitted files can be found
        # OutputPath = os.path.abspath(self.config['output']['output_path'])

        # Open the fitted h5 file
        with tables.open_file(h5name,'r+') as h5:
            node_path = '/ProcessingParams'
            h5.create_group(node_path,'ComputerInfo',title='Processing Computer Information',createparents=True)
            h5.create_group(node_path,'ConfigFiles',title='Config File Information',createparents=True)
            h5.create_array(node_path,'SoftwareVersion',np.array(Version),title='Version of software that made this file',createparents=True)
            h5.create_array(node_path,'RawFiles',np.array(RawFiles),title='The raw files used to produce this file',createparents=True)
            # h5.create_array(node_path,'OutputPath',np.array(OutputPath),title='Path where this file was originally made',createparents=True)
            node_path = '/ProcessingParams/ComputerInfo'
            h5.create_array(node_path,'PythonVersion',np.array(PythonVersion),title='Version of python used',createparents=True)
            h5.create_array(node_path,'Type',np.array(Type),title='Type of operating system',createparents=True)
            h5.create_array(node_path,'System',np.array(System),title='System information',createparents=True)
            h5.create_array(node_path,'User',np.array(User),title='Username',createparents=True)
            h5.create_array(node_path,'Host',np.array(Hostname),title='Hostname of the computer',createparents=True)
            node_path = '/ProcessingParams/ConfigFiles/File1'
            h5.create_array(node_path,'Name',np.array(Name),createparents=True)
            h5.create_array(node_path,'Path',np.array(Path),createparents=True)
            h5.create_array(node_path,'Contents',np.array(Contents),createparents=True)



# Don't think we need this if we save arrays as carrays and compress to start
# Class for repacking h5 files and compressing them
ignore_attribs = ['TITLE','CLASS','VERSION']
class repackh5(object):
    """
    Repack the input hdf5 file and compress the arrays
    """
    def __init__(self,input_fname,output_fname):
        self.input_file = input_fname
        self.output_file = output_fname

    def repack(self):
        FILTERS = tables.Filters(complib='zlib', complevel=1)
        with tables.open_file(self.input_file,'r') as input_h5:
            with tables.open_file(self.output_file,'w',filters=FILTERS) as output_h5:

                for group in input_h5.walk_groups():
                    group_name = group._v_pathname

                    # First make sure the group in the input file is in the output file.
                    output_h5_groups = [g._v_pathname for g in output_h5.walk_groups()]
                    if group_name not in output_h5_groups:
                        root, name = os.path.split(group_name)
                        output_h5.create_group(root,name)

                    # Now list the nodes in the group. For any node that isn't a group, write it
                    # to the output file.
                    for node in input_h5.list_nodes(group_name):
                        if node._v_attrs.CLASS != 'GROUP':
                            # Read the array, get it's attributes
                            name = node.name
                            title = node.title
                            data = node.read()
                            try:
                                shape = data.shape
                                is_array = True
                            except AttributeError:
                                shape = ()
                                is_array = False

                            #attributes
                            output_attribs = dict()
                            input_attribs = node.attrs.__dict__
                            output_attribs_keys = [x for x in input_attribs.keys() if not x.startswith('_') and not x in ignore_attribs]
                            for key in output_attribs_keys:
                                output_attribs[key] = input_attribs[key]


                            # Write the array, write it's attributes
                            if len(shape) and is_array:
                                atom = tables.Atom.from_dtype(data.dtype)
                                array = output_h5.create_carray(group_name,name,atom,shape,title=title)
                                array[:] = data
                            else:
                                output_h5.create_array(group_name,name,data,title=title)

                            for key in output_attribs_keys:
                                output_h5.set_node_attr(node._v_pathname,key,output_attribs[key])



def save_carray(h5, node, data, attributes):
    group, name = os.path.split(node)
    atom = tables.Atom.from_dtype(data.dtype)
    arr = h5.create_carray(group, name, atom, data.shape)
    arr[:] = data
    for k, v in attributes.items():
        h5.set_node_attr(node, k, v)



config_file_help = """Determine velocity vectors in altitude bins that are consistent with
input fitted line of sight velocities.

Requires a configuration file containing the following example format:
[DEFAULT]
[VVELS_OPTIONS]
#number of seconds of data to integrate
Recs2integrate=20
#use mean or median
mean_or_median=median
# altitude bin format: start,stop,stepsize,stride
# all in km, separate multiple bins settings with colon
altitude_bins=80,200,5,5:200,400,20,20
[INPUT]
# input paths (separate pulsetypes with comma,
# separate searches within same pulsetype with colons)
# example 2 pulsetypes, 2 search path per pulsetypes
file_paths=/path1/20150126.001_ac_5min-fitcal.h5:/path2/20150126.002_ac_5min-fitcal.h5,/path1/20150126.001_lp_5min-fitcal.h5:/path2/20150126.002_lp_5min-fitcal.h5
[OUTPUT]
# Output path
OUTPUT_PATH=/path/to/output/directory/%%(EXPNAME)s
# Output filename
OUTPUT_NAME=%%(OUTPUT_PATH)s/%%(EXPNAME)s_vvelsalt_%%(INTEG)s.h5
"""


# a function to run this file as a script
def main():
    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('config_file',help='A configuration file.')

    args = vars(parser.parse_args())
    vvelsalt = ResolveVectorsAlt(args['config_file'])
    vvelsalt.run()


# Run as a script
if __name__ == '__main__':
    main()
