import tables
from datetime import datetime
import numpy as np

class FittedVelocityDataHandler(object):
    def __init__(self,filelist):

        # Input checking
        # if not len(filelist):
        #     raise Exception('File list is empty.')
        self.filelist = filelist

        # # redundant error checking - consolidate
        # # check for FITS array in the files
        # has_fits = self.__has_fits_array(self.filelist[0])
        #
        # if not has_fits:
        #     raise Exception('Files are missing the /FittedParams/Fits array. Are you sure they are fitted data files?')

        # read times from the files
        # self.__catalog_data()

        # vvelsAlt does a lot of work checking arrays sizes
        # change to BTAFTP?
        # determine the array shapes in the files
        # self.__determine_array_sizes()

        # # initialize file handling variables
        # self.__loaded_file = None
        # self.__loaded_data = None
        # self.__loaded_beamcodes = None


    # def __has_fits_array(self,filename):
    #     with tables.open_file(filename,'r') as h5:
    #         try:
    #             node = h5.get_node('/FittedParams/Fits')
    #             return True
    #         except tables.NoSuchNodeError:
    #             return False


    # # This version has the ability to process multiple files together
    # # Should retain this and expand it to the vvelsLat
    # # Table this for later - kind of confusing how it works/should work
    # # For now, just manage one file at a time
    # def __catalog_data(self):
    #     # for each file in self.filelist, grab all the times and make a dictionary mapping to that file for each datetime
    #
    #     self.filetimes = dict()
    #
    #     # For every file that we have, grab the start and end times in the files.
    #     # We need to check every time, in case there is a gap in the file. Also
    #     # should make sure power, noise, cal, and unix time have same # of time records
    #     for fname in self.filelist:
    #         with tables.open_file(fname,'r') as h5:
    #             node = h5.get_node('/FittedParams')
    #             temp_times = h5.root.Time.UnixTime.read()
    #             fits_times = h5.root.FittedParams.Fits.shape[0]
    #             chi2_times = h5.root.FittedParams.FitInfo.chi2.shape[0]
    #             ne_times   = h5.root.FittedParams.Ne.shape[0]
    #             time_times = temp_times.shape[0]
    #
    #         num_times   = np.min([fits_times,chi2_times,ne_times,time_times])
    #
    #         for i in range(num_times):
    #             file_time = (datetime.utcfromtimestamp(temp_times[i,0]),
    #                          datetime.utcfromtimestamp(temp_times[i,1])
    #                         )
    #             self.filetimes[file_time] = fname
    #
    #     # now get an array of start and ends times from the keys of the filetimes dict
    #     temp = [list(x) for x in self.filetimes.keys()]
    #     temp.sort()
    #
    #     self.times = np.array(temp)
    #
    #
    # def __determine_array_sizes(self):
    #     temp_shape = dict()
    #     with tables.open_file(self.filelist[0],'r') as h5:
    #         temp_shape['chi2'] = h5.root.FittedParams.FitInfo.chi2.shape
    #         temp_shape['ne'] = h5.root.FittedParams.Ne.shape
    #         temp_shape['vel'] = h5.root.FittedParams.Fits.shape[0:3]
    #         temp_shape['altitude'] = h5.root.FittedParams.Altitude.shape
    #         temp_shape['mlat'] = h5.root.Geomag.MagneticLatitude.shape
    #
    #     self.__array_shapes = temp_shape

    def read_data(self, use_beams=None):

        self.data = {}

        # read data from standard AMISR fit files
        with tables.open_file(self.filelist,'r') as infile:

            # time
            self.data['time'] = infile.get_node('/Time/UnixTime')[:]

            # site
            lat = infile.get_node('/Site/Latitude').read()
            lon = infile.get_node('/Site/Longitude').read()
            alt = infile.get_node('/Site/Altitude').read()
            self.data['site'] = np.array([lat, lon, alt/1000.])

            # define which beams to use (default is all)
            self.data['beam_codes'] = infile.get_node('/BeamCodes')[:,0]
            if use_beams:
                bm_idx = np.array([i for i,b in enumerate(self.data['beam_codes']) if b in use_beams])
            else:
                bm_idx = np.arange(0,len(self.data['beam_codes']))

            # geodetic location of each measurement
            self.data['alt'] = infile.get_node('/Geomag/Altitude')[bm_idx,:].flatten()
            self.data['lat'] = infile.get_node('/Geomag/Latitude')[bm_idx,:].flatten()
            self.data['lon'] = infile.get_node('/Geomag/Longitude')[bm_idx,:].flatten()
            self.data['mlat'] = infile.get_node('/Geomag/MagneticLatitude')[bm_idx,:].flatten()

            # geodetic k vectors
            self.data['ke'] = infile.get_node('/Geomag/ke')[bm_idx,:].flatten()
            self.data['kn'] = infile.get_node('/Geomag/kn')[bm_idx,:].flatten()
            self.data['kz'] = infile.get_node('/Geomag/kz')[bm_idx,:].flatten()

            # ion masses
            ion_mass = infile.get_node('/FittedParams/IonMass')[:]
            nions = len(ion_mass)                                   # number of ions
            ion_idx = np.argwhere(ion_mass==16.).flatten()[0]       # find index for O+

            # line of sight velocity and error
            self.data['vlos'] = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,ion_idx,3].reshape((len(self.data['time'][:,0]),len(self.data['alt'])))
            self.data['dvlos'] = infile.get_node('/FittedParams/Errors')[:,bm_idx,:,ion_idx,3].reshape((len(self.data['time'][:,0]),len(self.data['alt'])))

            # chi2 and fitcode (for filtering poor quality data)
            self.data['chi2'] = infile.get_node('/FittedParams/FitInfo/chi2')[:,bm_idx,:].reshape((len(self.data['time'][:,0]),len(self.data['alt'])))
            self.data['fitcode'] = infile.get_node('/FittedParams/FitInfo/fitcode')[:,bm_idx,:].reshape((len(self.data['time'][:,0]),len(self.data['alt'])))

            # density (for filtering and ion upflow correction)
            self.data['ne'] = infile.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((len(self.data['time'][:,0]),len(self.data['alt'])))

            # # temperature (for ion upflow)
            # self.Te = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,-1,1].reshape((len(self.time[:,0]),len(self.alt)))
            # Ts = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:nions,1]
            # frac = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,:nions,0]
            # self.Ti = np.sum(Ts*frac,axis=-1).reshape((len(self.time[:,0]),len(self.alt)))
            #
            # # get up-B beam velocities for ion outflow correction
            # if self.upB_beamcode:
            #     upB_idx = np.argwhere(self.BeamCodes==self.upB_beamcode).flatten()
            #     if upB_idx:
            #         upB_alt = infile.get_node('/Geomag/Altitude')[upB_idx,:].flatten()
            #         upB_vlos = infile.get_node('/FittedParams/Fits')[:,upB_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(upB_alt)))
            #         upB_dvlos = infile.get_node('/FittedParams/Errors')[:,upB_idx,:,ion_idx,3].reshape((len(self.time[:,0]),len(upB_alt)))
            #         self.upB = {'alt':upB_alt, 'vlos':upB_vlos, 'dvlos':upB_dvlos}
            #     else:
            #         print('Warning: upB beam %d not found. Will not perform upflow subtraction.' % self.upB_beamcode)
            #         self.upB = None
            #         self.ionup = False

    def filter(self, ne=None, alt=None, mlat=None, chi2=None, fitcode=None, chirp=0.):

        with np.errstate(invalid='ignore'):

            # add chirp to LoS velocity
            self.data['vlos'] = self.data['vlos'] + chirp

            # discard data with low density
            if ne:
                inds = np.where((self.data['ne'] < ne[0]) | (self.data['ne'] > ne[1]))
                self.data['vlos'][inds] = np.nan
                self.data['dvlos'][inds] = np.nan

            # discard data outside of altitude range
            if alt:
                inds = np.where((self.data['alt'] < alt[0]*1000.) | (self.data['alt'] > alt[1]*1000.))
                self.data['vlos'][:,inds] = np.nan
                self.data['dvlos'][:,inds] = np.nan

            # discard data outside of magnetic latitude range
            if mlat:
                inds = np.where((self.data['mlat'] < mlat[0]) | (self.data['mlat'] > mlat[1]))
                self.data['vlos'][:,inds] = np.nan
                self.data['dvlos'][:,inds] = np.nan

            # discard data with extremely high or extremely low chi2 values
            if chi2:
                inds = np.where((self.data['chi2'] < chi2[0]) | (self.data['chi2'] > chi2[1]))
                self.data['vlos'][inds] = np.nan
                self.data['dvlos'][inds] = np.nan

            # discard data with poor fitcode (fitcodes 1-4 denote solution found, anything else should not be used)
            if fitcode:
                inds = np.where(~np.isin(self.data['fitcode'], fitcode))
                self.data['vlos'][inds] = np.nan
                self.data['dvlos'][inds] = np.nan


    # def __load_file(self,requested_time):
    #
    #     time_ind = np.where((self.times[:,0] <= requested_time) & (self.times[:,1] >= requested_time))[0]
    #     if len(time_ind) < 1:
    #         print("Requested time not found.")
    #         return None
    #
    #     # If it is, let's see if we have that file loaded or not
    #     needed_file = self.filetimes[tuple(self.times[time_ind[0]])]
    #     if self.__loaded_file != needed_file:
    #         print("Loading file: %s" % needed_file)
    #         # Load the needed arrays in FittedParams and the Time/UnixTime data
    #         with tables.open_file(needed_file,'r') as h5:
    #             # Handle beamcodes for AMISR and Sondrestrom
    #             try:
    #                 self.__loaded_beamcodes = np.array(h5.root.BeamCodes.read())
    #             except tables.NoSuchNodeError:
    #                 # SONDRESTROM NO SUPPORTED HERE. NEED TO ADD.
    #                 raise NotImplemented("Sondrestrom data not supported.")
    #                 # az = np.mean(output1['/Antenna']['Azimuth'])
    #                 # el = np.mean(output1['/Antenna']['Elevation'])
    #                 # self.__loaded_beamcodemap = np.array([[32768,az,el,0.0]])
    #
    #             # data
    #             temp_data = dict()
    #             temp_data['vel']  = h5.root.FittedParams.Fits[:,:,:,0,3]
    #             temp_data['evel'] = h5.root.FittedParams.Errors[:,:,:,0,3]
    #             temp_data['chi2'] = h5.root.FittedParams.FitInfo.chi2.read()
    #             temp_data['ne']   = h5.root.FittedParams.Ne.read()
    #             temp_data['mlat'] = h5.root.Geomag.MagneticLatitude.read()
    #
    #             # k-vectors
    #             temp_data['kvec'] = h5.root.Geomag.kvec.read()
    #
    #             # altitudes
    #             temp_data['altitude'] = h5.root.FittedParams.Altitude.read() / 1000.0
    #
    #             # site info
    #             temp_data['site'] = dict()
    #             temp_data['site']['lat'] = h5.root.Site.Latitude.read()
    #             temp_data['site']['lon'] = h5.root.Site.Longitude.read()
    #             temp_data['site']['alt'] = h5.root.Site.Altitude.read() / 1000.0
    #
    #             # Time
    #             temp_time = h5.root.Time.UnixTime.read()
    #
    #         self.__loaded_data = temp_data
        #
        #     # Check to make sure vel, chi2, ne, and unix_time
        #     # all have the same number of time records
        #     vel_times   = temp_data['vel'].shape[0]
        #     chi2_times  = temp_data['chi2'].shape[0]
        #     ne_times    = temp_data['ne'].shape[0]
        #     time_times  = temp_time.shape[0]
        #     num_times   = np.min([vel_times,chi2_times,ne_times,time_times])
        #
        #     # ensure all arrays only have the proper time array shape
        #     temp_data['vel']  = temp_data['vel'][:num_times,:]
        #     temp_data['evel'] = temp_data['evel'][:num_times,:]
        #     temp_data['chi2'] = temp_data['chi2'][:num_times,:]
        #     temp_data['ne']   = temp_data['ne'][:num_times,:]
        #     temp_time = temp_time[:num_times,:]
        #
        #     # Convert UnixTime to datetime
        #     temp = list()
        #     for i in range(num_times):
        #         time_pair = [datetime.utcfromtimestamp(temp_time[i,0]),
        #                      datetime.utcfromtimestamp(temp_time[i,1])]
        #         temp.append(time_pair)
        #     self.__loaded_time = np.array(temp)
        #     self.__loaded_file = needed_file
        #
        # return True


    # def get_record(self,requested_time):
    #     if self.__load_file(requested_time) is None:
    #         return None
    #
    #     time_ind = np.where((self.__loaded_time[:,0] <= requested_time) & (self.__loaded_time[:,1] >= requested_time))[0]
    #
    #     data = dict()
    #     data['vel']  = np.squeeze(self.__loaded_data['vel'][time_ind,:])
    #     data['evel'] = np.squeeze(self.__loaded_data['evel'][time_ind,:])
    #     data['chi2'] = np.squeeze(self.__loaded_data['chi2'][time_ind,:])
    #     data['ne']   = np.squeeze(self.__loaded_data['ne'][time_ind,:])
    #
    #     data['kvec']     = self.__loaded_data['kvec']
    #     data['altitude'] = self.__loaded_data['altitude']
    #     data['mlat'] = self.__loaded_data['mlat']
    #
    #     data['site'] = self.__loaded_data['site']
    #
    #     data['bmcodes'] = self.__loaded_beamcodes
    #
    #     return data

    def get_records(self,start_time,end_time):

        # figure out how many time records we have to get
        # the logic on the next line is correct, even though it seems confusing at first
        # The start time needs to be checked against the end times of each time record
        # and the e time needs to be checked against the start times of each record.
        request_time_inds = np.where((self.data['time'][:,0] < end_time) & (self.data['time'][:,1] >= start_time))[0]
        temp_times = self.data['time'][request_time_inds]
        # print(request_time_inds, temp_times.shape, self.data['vlos'].shape, self.data['vlos'][request_time_inds].shape)

        data = {'vel':self.data['vlos'][request_time_inds],
                'evel':self.data['dvlos'][request_time_inds],
                'chi2':self.data['chi2'][request_time_inds],
                'ne':self.data['ne'][request_time_inds]}

        # if len(temp_times) < 1:
        #     print("No data for request start and end times.")
        #     return None
        #
        # request_times = list()
        # epoch = datetime(1970,1,1)
        # for tup in temp_times:
        #     t1 = (tup[0] - epoch).total_seconds()
        #     t2 = (tup[1] - epoch).total_seconds()
        #     request_times.append(datetime.utcfromtimestamp((t1 + t2) / 2.0))
        # request_times.sort()
        #
        # num_times = len(request_times)
        #
        # arrsh = self.__array_shapes
        # data = dict()
        # data['vel']  = np.zeros((num_times,) + arrsh['vel'][1:])
        # data['evel'] = np.zeros((num_times,) + arrsh['vel'][1:])
        # data['chi2'] = np.zeros((num_times,) + arrsh['chi2'][1:])
        # data['ne']   = np.zeros((num_times,) + arrsh['ne'][1:])
        #
        # # now get the data for the requested time
        # for i,time in enumerate(request_times):
        #     temp = self.get_record(time)
        #     if i == 0:
        #         data['kvec']     = temp['kvec']
        #         data['altitude'] = temp['altitude']
        #         data['site']     = temp['site']
        #         data['bmcodes']  = temp['bmcodes']
        #         data['mlat']     = temp['mlat']
        #
        #     data['vel'][i,:]  = temp['vel']
        #     data['evel'][i,:] = temp['evel']
        #     data['chi2'][i,:] = temp['chi2']
        #     data['ne'][i,:]   = temp['ne']

        return data, temp_times
