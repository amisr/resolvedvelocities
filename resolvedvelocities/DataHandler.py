import tables
from datetime import datetime
import numpy as np

class FittedVelocityDataHandler(object):
    def __init__(self,filename):

        self.filename = filename



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


    def load_data(self, use_beams=None):

        self.data = {}

        # read data from standard AMISR fit files
        with tables.open_file(self.filename, 'r') as infile:

            # time
            self.utime = infile.get_node('/Time/UnixTime')[:]
            self.Nrecords = len(self.utime[:,0])

            # site
            lat = infile.get_node('/Site/Latitude').read()
            lon = infile.get_node('/Site/Longitude').read()
            alt = infile.get_node('/Site/Altitude').read()
            self.site = np.array([lat, lon, alt/1000.])

            # define which beams to use (default is all)
            self.beam_codes = infile.get_node('/BeamCodes')[:,0]
            if use_beams:
                bm_idx = np.array([i for i,b in enumerate(self.beam_codes) if b in use_beams])
            else:
                bm_idx = np.arange(0,len(self.beam_codes))

            # geodetic location of each measurement
            self.alt = infile.get_node('/Geomag/Altitude')[bm_idx,:].flatten()/1000.
            self.lat = infile.get_node('/Geomag/Latitude')[bm_idx,:].flatten()
            self.lon = infile.get_node('/Geomag/Longitude')[bm_idx,:].flatten()
            self.mlat = infile.get_node('/Geomag/MagneticLatitude')[bm_idx,:].flatten()

            self.Nbins = len(self.alt)

            # geodetic k vectors
            self.ke = infile.get_node('/Geomag/ke')[bm_idx,:].flatten()
            self.kn = infile.get_node('/Geomag/kn')[bm_idx,:].flatten()
            self.kz = infile.get_node('/Geomag/kz')[bm_idx,:].flatten()

            # ion masses
            ion_mass = infile.get_node('/FittedParams/IonMass')[:]
            nions = len(ion_mass)                                   # number of ions
            ion_idx = np.argwhere(ion_mass==16.).flatten()[0]       # find index for O+

            # line of sight velocity and error
            self.vlos = infile.get_node('/FittedParams/Fits')[:,bm_idx,:,ion_idx,3].reshape((self.Nrecords, self.Nbins))
            self.dvlos = infile.get_node('/FittedParams/Errors')[:,bm_idx,:,ion_idx,3].reshape((self.Nrecords, self.Nbins))

            # chi2 and fitcode (for filtering poor quality data)
            self.chi2 = infile.get_node('/FittedParams/FitInfo/chi2')[:,bm_idx,:].reshape((self.Nrecords, self.Nbins))
            self.fitcode = infile.get_node('/FittedParams/FitInfo/fitcode')[:,bm_idx,:].reshape((self.Nrecords, self.Nbins))

            # density (for filtering and ion upflow correction)
            self.ne = infile.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((self.Nrecords, self.Nbins))

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
            self.vlos = self.vlos + chirp

            # discard data with low density
            if ne:
                inds = np.where((self.ne < ne[0]) | (self.ne > ne[1]))
                self.vlos[inds] = np.nan
                self.dvlos[inds] = np.nan

            # discard data outside of altitude range
            if alt:
                inds = np.where((self.alt < alt[0]) | (self.alt > alt[1]))
                self.vlos[:,inds] = np.nan
                self.dvlos[:,inds] = np.nan

            # discard data outside of magnetic latitude range
            if mlat:
                inds = np.where((self.mlat < mlat[0]) | (self.mlat > mlat[1]))
                self.vlos[:,inds] = np.nan
                self.dvlos[:,inds] = np.nan

            # discard data with extremely high or extremely low chi2 values
            if chi2:
                inds = np.where((self.chi2 < chi2[0]) | (self.chi2 > chi2[1]))
                self.vlos[inds] = np.nan
                self.dvlos[inds] = np.nan

            # discard data with poor fitcode (fitcodes 1-4 denote solution found, anything else should not be used)
            if fitcode:
                inds = np.where(~np.isin(self.fitcode, fitcode))
                self.vlos[inds] = np.nan
                self.dvlos[inds] = np.nan


    def get_record_indices(self, start_time, end_time):
        request_time_inds = np.where((self.utime[:,0] < end_time) & (self.utime[:,1] > start_time))[0]
        return request_time_inds


    # def get_records(self,start_time,end_time):
    #
    #     # figure out how many time records we have to get
    #     # the logic on the next line is correct, even though it seems confusing at first
    #     # The start time needs to be checked against the end times of each time record
    #     # and the e time needs to be checked against the start times of each record.
    #     request_time_inds = np.where((self.utime[:,0] < end_time) & (self.utime[:,1] > start_time))[0]
    #     temp_times = self.data['time'][request_time_inds]
    #     # print(request_time_inds, temp_times.shape, self.data['vlos'].shape, self.data['vlos'][request_time_inds].shape)
    #
    #     data = {'vel':self.data['vlos'][request_time_inds],
    #             'evel':self.data['dvlos'][request_time_inds],
    #             'chi2':self.data['chi2'][request_time_inds],
    #             'ne':self.data['ne'][request_time_inds]}

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
