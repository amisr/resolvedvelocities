import numpy as np
import datetime as dt
from apexpy import Apex
from scipy import interpolate
import coord_convert as cc
import tables
import sys
import os


class Field(object):
    def __init__(self, lat, lon, field, alt):

        self.map_velocity_field(lat.flatten(), lon.flatten(), field.reshape(-1,field.shape[-1]), alt)
        self.convert_to_ECEF()
        self.create_interpolators()


    def map_velocity_field(self, alat, alon, field, alt_in):
        # alat - 1D array (N,) of apex magnetic latitude
        # alon - 1D array (N,) of apex magnetic longitude
        # field - array (N,3) of geodetic E, N, U components of the velocity field at each position
        # alt_in - scalar geodetic altitude of specified points

        # define output altitudes
        altitude = np.arange(50., 1000., 50.)

        # initialize Apex object
        A = Apex(date=2019)

        # find positions at each altitude
        self.altitude = np.repeat(altitude,len(alat))
        self.latitude, self.longitude, __ = A.apex2geo(np.tile(alat,len(altitude)), np.tile(alon,len(altitude)), self.altitude)

        # map field to each altitude
        f = np.array([A.map_V_to_height(alat, alon, alt_in, a, field.T).T for a in altitude])
        self.field = f.reshape(-1,f.shape[-1])


    def convert_to_ECEF(self):

        self.X, self.Y, self.Z = cc.geodetic_to_cartesian(self.latitude, self.longitude, self.altitude)
        self.Vx, self.Vy, self.Vz = cc.vector_geodetic_to_cartesian(self.field[:,1], self.field[:,0], self.field[:,2], self.latitude, self.longitude, self.altitude)

    def create_interpolators(self):

        self.interpVx = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vx)
        self.interpVy = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vy)
        self.interpVz = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vz)


class Radar(object):
    def __init__(self, site, beams=None, azimuth=None, elevation=None, range_step=50.):
        # beam - list of beam codes
        # azimuth - list of azimuth angle for each beam (in degrees)
        # elevation - list of elevation angle for each beam (in degrees)
        # range_step - step between each range gate (in km)

        self.site = site
        if beams:
            bc = np.loadtxt(os.path.join(os.path.dirname(__file__), 'bcotable.txt'))
            idx = np.where(np.in1d(bc[:,0],beams))[0]
            self.beam_codes = bc[idx,:]
        elif azimuth:
            self.beam_codes = np.array([range(len(azimuth)),azimuth,elevation,np.full(len(azimuth),np.nan)]).T

        self.range_step = range_step
        self.X0, self.Y0, self.Z0 = cc.geodetic_to_cartesian(site[0], site[1], site[2])
        self.get_gate_locations(self.beam_codes[:,1], self.beam_codes[:,2], range_step)
        self.geodetic_locations()

    def get_gate_locations(self, az, el, rs):

        # create array of ranges
        ranges = np.arange(80.,800., rs)*1000.

        # find E, N, U components of k (find k vector in geodetic coordinates)
        az = np.array(az)*np.pi/180.
        el = np.array(el)*np.pi/180.
        ke = np.cos(el)*np.sin(az)
        kn = np.cos(el)*np.cos(az)
        ku = np.sin(el)

        # convert geodetic k vector to ECEF
        kx, ky, kz = cc.vector_geodetic_to_cartesian(kn, ke, ku, self.site[0], self.site[1], self.site[2])

        # calculate position of each range gate in ECEF
        self.X = kx[:,None]*ranges + self.X0
        self.Y = ky[:,None]*ranges + self.Y0
        self.Z = kz[:,None]*ranges + self.Z0

        # form array of k vector at each range gate
        self.kvec = np.array([np.tile(np.array([x,y,z]), (len(ranges),1)) for x, y, z in zip(kx,ky,kz)])


    def geodetic_locations(self):
        # calculate gate position and k vectors in geodetic coordinates
        self.lat, self.lon, self.alt = cc.cartesian_to_geodetic(self.X, self.Y, self.Z)
        self.kn, self.ke, self.kz = cc.vector_cartesian_to_geodetic(self.kvec[:,:,0], self.kvec[:,:,1], self.kvec[:,:,2], self.X, self.Y, self.Z)



def create_dataset(field, radar):
    # input: Field object, Radar object

    # interpolate the field to the radar bin locations
    Vx = field.interpVx(np.array([radar.X, radar.Y, radar.Z]).T)
    Vy = field.interpVy(np.array([radar.X, radar.Y, radar.Z]).T)
    Vz = field.interpVz(np.array([radar.X, radar.Y, radar.Z]).T)
    Vvec = np.array([Vx, Vy, Vz]).T

    # create unix time array
    time0 = (dt.datetime(2019,5,28,0,0)-dt.datetime.utcfromtimestamp(0)).total_seconds()
    times = np.array([[time0+t*60., time0+(t+1)*60.] for t in range(10)])

    # calculate LoS velocity for each bin by taking the dot product of the radar kvector and the interpolated field
    Vlos = np.tile(np.einsum('...i,...i->...',radar.kvec, Vvec), (len(times),1,1))
    # assume constant error
    dVlos = np.full(Vlos.shape, 10.)

    # create fit and error arrays that match the shape of whats in the processed fitted files
    s = Vlos.shape
    fit_array = np.full((s[0],s[1],s[2],6,4),np.nan)
    fit_array[:,:,:,0,3] = Vlos
    err_array = np.full((s[0],s[1],s[2],6,4),np.nan)
    err_array[:,:,:,0,3] = dVlos

    # generate dummy density array
    ne = np.full(Vlos.shape, 1e11)

    # create output hdf5 file
    filename = 'synthetic_data.h5'
    with tables.open_file(filename, mode='w') as file:
        file.create_group('/','FittedParams')
        file.create_group('/','Geomag')
        file.create_group('/','Time')
        file.create_group('/','Site')

        file.create_array('/','BeamCodes',radar.beam_codes)

        file.create_array('/FittedParams','Fits',fit_array)
        file.create_array('/FittedParams','Errors',err_array)
        file.create_array('/FittedParams','Ne',ne)

        file.create_array('/Geomag','Latitude',radar.lat)
        file.create_array('/Geomag','Longitude',radar.lon)
        file.create_array('/Geomag','Altitude',radar.alt*1000.)

        file.create_array('/Geomag','ke',radar.ke)
        file.create_array('/Geomag','kn',radar.kn)
        file.create_array('/Geomag','kz',radar.kz)

        file.create_array('/Time','UnixTime',times)

        file.create_array('/Site','Latitude',radar.site[0])
        file.create_array('/Site','Longitude',radar.site[1])
        file.create_array('/Site','Altitude',radar.site[2])