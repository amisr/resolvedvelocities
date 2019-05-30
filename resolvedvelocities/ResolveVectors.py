# ResolveVectors.py

import tables
import numpy as np
import datetime as dt
import configparser
from apexpy import Apex


class ResolveVectors(object):
    def __init__(self, config=None):
        self.read_config(config)

    def read_config(self, config_file):
        # read config file
        config = configparser.ConfigParser(allow_no_values=True)
        config.read(config_file)

        self.datafile = config.get('DEFAULT', 'DATAFILE')
        self.chirp = eval(config.get('DEFAULT', 'CHIRP'))
        self.neMin = eval(config.get('DEFAULT', 'NEMIN'))
        self.integration_time = config.getfloat('DEFAULT', 'INTTIME', fallback=None)
        self.covar = eval(config.get('DEFAULT', 'COVAR'))
        self.ppp = eval(config.get('DEFAULT', 'PPP'))
        self.minalt = eval(config.get('DEFAULT', 'MINALT'))
        self.maxalt = eval(config.get('DEFAULT', 'MAXALT'))
        self.minnumpoints = eval(config.get('DEFAULT', 'MINNUMPOINTS'))
        self.upB_beamcode = config.getint('DEFAULT', 'UPB_BEAMCODE', fallback=None)
        self.ionup = config.get('DEFAULT', 'IONUP', fallback=None)
        self.use_beams = config.get('DEFAULT', 'USE_BEAMS', fallback=None)
        if self.use_beams:
            self.use_beams = eval(self.use_beams)

        # list of beam codes to use

    def read_data(self):
        # read data from standard AMISR fit files
        with tables.open_file(self.datafile,'r') as file:

            # time
            self.time = file.get_node('/Time/UnixTime')[:]

            # beam codes
            # define which beams to use (default is all)
            self.BeamCodes=file.get_node('/BeamCodes')[:,0]
            if self.use_beams:
                bm_idx = np.array([i for i,b in enumerate(self.BeamCodes) if b in self.use_beams])
            else:
                bm_idx = np.arange(0,len(self.BeamCodes))
            print bm_idx

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
            self.dvlos = file.get_node('/FittedParams/Errors')[:,bm_idx,:,0,3].reshape((len(self.time[:,0]),len(self.alt)))

            # density (for filtering and ion upflow)
            self.ne = file.get_node('/FittedParams/Ne')[:,bm_idx,:].reshape((len(self.time[:,0]),len(self.alt)))

            # temperature (for ion upflow)
            self.Te = file.get_node('/FittedParams/Fits')[:,bm_idx,:,5,1].reshape((len(self.time[:,0]),len(self.alt)))
            Ts = file.get_node('/FittedParams/Fits')[:,bm_idx,:,:5,1]
            frac = file.get_node('/FittedParams/Fits')[:,bm_idx,:,:5,0]
            self.Ti = np.sum(Ts*frac,axis=-1).reshape((len(self.time[:,0]),len(self.alt)))

            # get up-B beam velocities for ion outflow correction
            if self.upB_beamcode:
                upB_idx = np.argwhere(self.BeamCodes==self.upB_beamcode).flatten()
                upB_alt = file.get_node('/Geomag/Altitude')[upB_idx,:].flatten()
                upB_vlos = file.get_node('/FittedParams/Fits')[:,upB_idx,:,0,3].reshape((len(self.time[:,0]),len(upB_alt)))
                upB_dvlos = file.get_node('/FittedParams/Errors')[:,upB_idx,:,0,3].reshape((len(self.time[:,0]),len(upB_alt)))
                self.upB = {'alt':upB_alt, 'vlos':upB_vlos, 'dvlos':upB_dvlos}



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

        # find indices where nans will be removed and should be inserted in new arrays
        replace_nans = np.array([r-i for i,r in enumerate(np.argwhere(np.isnan(self.alt)).flatten())])

        glat = self.lat[np.isfinite(self.lat)]
        glon = self.lon[np.isfinite(self.lon)]
        galt = self.alt[np.isfinite(self.alt)]/1000.

        # intialize apex coordinates
        self.Apex = Apex(date=dt.datetime.utcfromtimestamp(self.time[0,0]))

        # find magnetic latitude and longitude
        mlat, mlon = self.Apex.geo2apex(glat, glon, galt)
        self.mlat = np.insert(mlat,replace_nans,np.nan)
        self.mlon = np.insert(mlon,replace_nans,np.nan)

        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(glat, glon, galt)
        d1 = np.insert(d1,replace_nans,np.nan,axis=1)
        d2 = np.insert(d2,replace_nans,np.nan,axis=1)
        d3 = np.insert(d3,replace_nans,np.nan,axis=1)
        d = np.array([d1,d2,d3]).T

        # kvec in geodetic coordinates [e n u]
        kvec = np.array([self.ke, self.kn, self.kz]).T

        # find components of k for e1, e2, e3 basis vectors (Laundal and Richmond, 2016 eqn. 60)
        self.A = np.einsum('ij,ijk->ik', kvec, d)

        # calculate scaling factor D, used for ion outflow correction (Richmond, 1995 eqn. 3.15)
        d1_cross_d2 = np.cross(d1.T,d2.T).T
        self.D = np.sqrt(np.sum(d1_cross_d2**2,axis=0))



    def ion_upflow_correction(self):

        # correct the los velocities for the entire array at each time
        for t in range(len(self.time)):
            if self.ionup == 'UPB':
                # interpolate velocities from up B beam to all other measurements 
                vion, dvion = lin_interp(self.alt, self.upB['alt'], self.upB['vlos'][t], self.upB['dvlos'][t])
            elif self.ionup == 'EMP':
                # use empirical method to find ion upflow
                # NOTE: NOT DEVELOPED YET!!!
                vion, dvion = ion_upflow(self.Te, self.Ti, self.ne, self.alt)

            # LoS velocity correction to remove ion upflow
            self.vlos[t] = self.vlos[t] + self.D*self.A[:,2]*vion
            # corrected error in new LoS velocities
            self.dvlos[t] = np.sqrt(self.dvlos[t]**2 + self.D**2*self.A[:,2]**2*dvion**2)



    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in some way by initial config file
        # each bin has a specified MLAT/MLON

        bin_edge_mlat = np.arange(64.0,68.0,0.25)
        self.bin_mlat = (bin_edge_mlat[:-1]+bin_edge_mlat[1:])/2.
        self.bin_mlon = np.full(self.bin_mlat.shape, np.nanmean(self.mlon))
        self.bin_idx = [np.argwhere((self.mlat>=bin_edge_mlat[i]) & (self.mlat<bin_edge_mlat[i+1])).flatten() for i in range(len(bin_edge_mlat)-1)]



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


    def compute_vectors(self):
        # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
        
        Velocity = []
        VelocityCovariance = []

        # For each integration period and bin, calculate covarient components of drift velocity (Ve1, Ve2, Ve3)
        # loop over integration periods
        for tidx in self.int_idx:
            Vel = []
            SigmaV = []
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
        self.bin_glat, self.bin_glon, err = self.Apex.apex2geo(self.bin_mlat, self.bin_mlon, alt)
        # self.gdlat = gdlat
        # self.gdlon = gdlon
        self.bin_galt = np.full(self.bin_glat.shape, alt)

        # apex basis vectors in geodetic coordinates [e n u]
        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.Apex.basevectors_apex(self.bin_mlat, self.bin_mlon, alt, coords='apex')

        e = np.array([e1,e2,e3]).T
        self.Velocity_gd = np.einsum('ijk,...ik->...ij',e,self.Velocity)
        self.VelocityCovariance_gd = np.einsum('ijk,...ikl,iml->...ijm',e,self.VelocityCovariance,e)

        d = np.array([d1,d2,d3]).T
        self.ElectricField_gd = np.einsum('ijk,...ik->...ij',d,self.ElectricField)
        self.ElectricFieldCovariance_gd = np.einsum('ijk,...ikl,iml->...ijm',d,self.ElectricFieldCovariance,d)

        # calculate vector magnitude and direction
        self.Vgd_mag, self.Vgd_mag_err, self.Vgd_dir, self.Vgd_dir_err = magnitude_direction(self.Velocity_gd, self.VelocityCovariance_gd, -e2.T)
        self.Egd_mag, self.Egd_mag_err, self.Egd_dir, self.Egd_dir_err = magnitude_direction(self.ElectricField_gd, self.ElectricFieldCovariance_gd, -e2.T)


    def save_output(self):
        # TODO: come up with a better way to manage all this
        # - naming conventions?  Do these make sense?
        # - add processing paramters
        # - add config file info
        # - add computer/user info
        
        # save output file
        filename = 'test_vvels.h5'

        with tables.open_file(filename,mode='w') as file:

            # copy some groups directly from fitted input file
            with tables.open_file(self.datafile, mode='r') as inputfile:
                if not self.integration_time:
                    file.copy_children(inputfile.get_node('/Time'), file.create_group('/','Time'))
                file.copy_children(inputfile.get_node('/Site'), file.create_group('/','Site'))

            file.create_group('/', 'Magnetic')

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

            file.create_group('/', 'Geographic')

            file.create_array('/Geographic', 'GeographicLatitude', self.bin_glat)
            file.set_node_attr('/Geographic/GeographicLatitude', 'TITLE', 'Geographic Latitude')
            file.set_node_attr('/Geographic/GeographicLatitude', 'Size', 'Nbins')

            file.create_array('/Geographic','GeographicLongitude', self.bin_glon)
            file.set_node_attr('/Geographic/GeographicLongitude', 'TITLE', 'Geographic Longitude')
            file.set_node_attr('/Geographic/GeographicLongitude', 'Size', 'Nbins')

            file.create_array('/Geographic', 'Velocity', self.Velocity_gd)
            file.set_node_attr('/Geographic/Velocity', 'TITLE', 'Plama Drift Velocity')
            file.set_node_attr('/Geographic/Velocity', 'Size', 'Nrecords x Nbins x 3 (East, North, Up)')
            file.set_node_attr('/Geographic/Velocity', 'Units', 'm/s')

            file.create_array('/Geographic','SigmaV', self.VelocityCovariance_gd)
            file.set_node_attr('/Geographic/SigmaV', 'TITLE', 'Velocity Covariance Matrix')
            file.set_node_attr('/Geographic/SigmaV', 'Size', 'Nrecords x Nbins x 3 x 3')
            file.set_node_attr('/Geographic/SigmaV', 'Units', 'm/s')

            file.create_array('/Geographic','Vmag',self.Vgd_mag)
            file.set_node_attr('/Geographic/Vmag', 'TITLE', 'Velocity Magnitude')
            file.set_node_attr('/Geographic/Vmag', 'Size', 'Nrecords x Nbins')
            file.set_node_attr('/Geographic/Vmag', 'Units', 'm/s')

            file.create_array('/Geographic','errVmag',self.Vgd_mag_err)
            file.set_node_attr('/Geographic/errVmag', 'TITLE', 'Velocity Magnitude Error')
            file.set_node_attr('/Geographic/errVmag', 'Size', 'Nrecords x Nbins')
            file.set_node_attr('/Geographic/errVmag', 'Units', 'm/s')

            file.create_array('/Geographic','Vdir',self.Vgd_dir)
            file.set_node_attr('/Geographic/Vdir', 'TITLE', 'Velocity Direction from Magnetic North Meridian')
            file.set_node_attr('/Geographic/Vdir', 'Size', 'Nrecord x Nbins')
            file.set_node_attr('/Geographic/Vdir', 'Units', 'Degrees')

            file.create_array('/Geographic','errVdir',self.Vgd_dir_err)
            file.set_node_attr('/Geographic/errVdir', 'TITLE', 'Error in Velocity Direction from Magnetic North Meridian')
            file.set_node_attr('/Geographic/errVdir', 'Size', 'Nrecord x Nbins')
            file.set_node_attr('/Geographic/errVdir', 'Units', 'Degrees')

            file.create_array('/Geographic','ElectricField',self.ElectricField_gd)
            file.set_node_attr('/Geographic/ElectricField', 'TITLE', 'Convection Electric Field')
            file.set_node_attr('/Geographic/ElectricField', 'Size', 'Nrecords x Nbins x 3 (East, North, Up)')
            file.set_node_attr('/Geographic/ElectricField', 'Units', 'V/m')

            file.create_array('/Geographic','SigmaE',self.ElectricFieldCovariance_gd)
            file.set_node_attr('/Geographic/SigmaE', 'TITLE', 'Electric Field Covariance Matrix')
            file.set_node_attr('/Geographic/SigmaE', 'Size', 'Nrecords x Nbins x 3 x 3')
            file.set_node_attr('/Geographic/SigmaE', 'Units', 'V/m')

            file.create_array('/Geographic','Emag',self.Egd_mag)
            file.set_node_attr('/Geographic/Emag', 'TITLE', 'Electric Field Magnitude')
            file.set_node_attr('/Geographic/Emag', 'Size', 'Nrecords x Nbins')
            file.set_node_attr('/Geographic/Emag', 'Units', 'V/m')

            file.create_array('/Geographic','errEmag',self.Egd_mag_err)
            file.set_node_attr('/Geographic/errEmag', 'TITLE', 'Electric Field Magnitude Error')
            file.set_node_attr('/Geographic/errEmag', 'Size', 'Nrecords x Nbins')
            file.set_node_attr('/Geographic/errEmag', 'Units', 'V/m')

            file.create_array('/Geographic','Edir',self.Egd_dir)
            file.set_node_attr('/Geographic/Edir', 'TITLE', 'Electric Field Direction from Magnetic North Meridian')
            file.set_node_attr('/Geographic/Edir', 'Size', 'Nrecord x Nbins')
            file.set_node_attr('/Geographic/Edir', 'Units', 'Degrees')

            file.create_array('/Geographic','errEdir',self.Egd_dir_err)
            file.set_node_attr('/Geographic/errEdir', 'TITLE', 'Error in Electric Field Direction from Magnetic North Meridian')
            file.set_node_attr('/Geographic/errEdir', 'Size', 'Nrecord x Nbins')
            file.set_node_attr('/Geographic/errEdir', 'Units', 'Degrees')



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


def lin_interp(x, xp, fp, dfp):
    # Piecewise linear interpolation routine that returns interpolated values and their errors

    # find the indicies of xp that bound each value in x
    # Note: where x is out of range of xp, -1 is used as a place holder
    #   This provides a valid "dummy" index for the array calculations and can be used to identify values to nan in final output
    i = np.array([np.argwhere((xi>=xp[:-1]) & (xi<xp[1:])).flatten()[0] if ((xi>=np.nanmin(xp)) & (xi<np.nanmax(xp))) else -1 for xi in x])
    # calculate X
    X = (x-xp[i])/(xp[i+1]-xp[i])
    # calculate interpolated values
    f = fp[i] + (fp[i+1]-fp[i])*X
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
    # This is all done with somewhat obtuse matrix algebra using einsum to prevent nested for loops
    # Input vectors are assumed to have orthogonal components

    AA = np.einsum('...i,...i->...', A, A)                  # dot product of A and A
    ASA = np.einsum('...i,...ij,...j->...', A, Sig, A)      # matrix multipy A*Sig*A
    ee = np.einsum('...i,...i->...', e, e)                  # dot product of e and e
    eA = np.einsum('...i,...i->...', e, A)                  # dot product of e and A

    magnitude = np.sqrt(AA)
    mag_err = np.sqrt(ASA/AA)

    # find ep, perpendicular to both e and geodetic up
    ep = np.cross(e,np.array([0,0,1]))
    epep = np.einsum('...i,...i->...', ep, ep)
    epA = np.einsum('...i,...i->...', ep, A)

    direction = np.arctan2(np.sqrt(ee)*epA,np.sqrt(epep)*eA)

    B = np.einsum('ij,...i->...ij',ep,eA)-np.einsum('ij,...i->...ij',e,epA)     # B = ep(e*A)-e(ep*A) = A x (ep x e)
    BSB = np.einsum('...i,...ij,...j->...', B, Sig, B)      # matrix multipy B*Sig*B
    dir_err = np.sqrt(epep*ee*BSB)/(ee*epA**2-epep*eA**2)

    return magnitude, mag_err, direction*180./np.pi, dir_err*180./np.pi

