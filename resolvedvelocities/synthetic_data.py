import numpy as np
import datetime as dt
from apexpy import Apex
from scipy import interpolate
import coord_convert as cc
import tables


class Field(object):
    def __init__(self, lat, lon, field, alt):

        self.map_velocity_field([lat, lon], field, alt)

        self.convert_to_ECEF()
        self.create_interpolators()


    def map_velocity_field(self, grid, field, alt_in):
        # input: apex grid at altitude alt_in
        # geodetic components of velocity on grid
        # output: 3D geodetic latitude, longitude, altitude
        # 3D grid of geodetic velocity components

        # define output altitudes
        altitude = np.arange(50., 1000., 50.)

        A = Apex(date=2019)

        field_out = []
        lat_out = []
        lon_out = []
        alt_out = []
        fs = field.shape
        for a in altitude:

            full_field = A.map_V_to_height(grid[0].flatten(), grid[1].flatten(), alt_in, a, field.reshape((fs[0]*fs[1],fs[2])).T)
            field_out.append(full_field.T)

            lat, lon, __ = A.apex2geo(grid[0].flatten(),grid[1].flatten(),a)
            lat_out.append(lat)
            lon_out.append(lon)
            alt_out.append(np.full(lat.shape,a))

        self.field = np.array(field_out).reshape((fs[0]*fs[1]*len(altitude),fs[2]))
        self.latitude = np.array(lat_out).flatten()
        self.longitude = np.array(lon_out).flatten()
        self.altitude = np.array(alt_out).flatten()

    def convert_to_ECEF(self):

        self.X, self.Y, self.Z = cc.geodetic_to_cartesian(self.latitude, self.longitude, self.altitude)
        self.Vx, self.Vy, self.Vz = cc.vector_geodetic_to_cartesian(self.field[:,1], self.field[:,0], self.field[:,2], self.latitude, self.longitude, self.altitude)

    def create_interpolators(self):

        self.interpVx = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vx)
        self.interpVy = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vy)
        self.interpVz = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vz)


class Radar(object):
    def __init__(self, site, azimuth, elevation, range_step):

        self.site = site
        self.azimuth = azimuth
        self.elevation = elevation
        self.range_step = range_step
        self.X0, self.Y0, self.Z0 = cc.geodetic_to_cartesian(site[0], site[1], site[2])
        self.bin_locations(azimuth, elevation, range_step)
        self.geodetic_locations()

    def bin_locations(self, az, el, rs):
        # a lot of this can probably be handled more elegantly with clever broadcasting, but this works for now

        ranges = np.arange(80.,800., rs)*1000.

        self.X = []
        self.Y = []
        self.Z = []
        self.kvec = []
        self.ke = []
        self.kn = []
        self.kz = []

        for a, e in zip(az,el):
            ke = np.cos(e*np.pi/180.)*np.sin(a*np.pi/180.)
            kn = np.cos(e*np.pi/180.)*np.cos(a*np.pi/180.)
            ku = np.sin(e*np.pi/180.)
            # self.ke.append(ke)
            # self.kn.append(kn)
            # self.kz.append(kz)

            kx, ky, kz = cc.vector_geodetic_to_cartesian(kn, ke, ku, self.site[0], self.site[1], self.site[2])

            X = kx*ranges + self.X0
            Y = ky*ranges + self.Y0
            Z = kz*ranges + self.Z0
            self.X.append(X)
            self.Y.append(Y)
            self.Z.append(Z)
            self.kvec.append(np.tile(np.array([kx, ky, kz]), (len(ranges),1)))

            # kn, ke, ku = cc.vector_cartesian_to_geodetic(kx, ky, kz, X, Y, Z)
            # self.ke.append(ke)
            # self.kn.append(kn)
            # self.kz.append(kz)

        self.X = np.array(self.X)
        self.Y = np.array(self.Y)
        self.Z = np.array(self.Z)
        self.kvec = np.array(self.kvec)
        # self.ke = np.array(self.ke)
        # self.kn = np.array(self.kn)
        # self.kz = np.array(self.kz)

        # print self.X.shape, self.kvec.shape


    def geodetic_locations(self):

        self.lat, self.lon, self.alt = cc.cartesian_to_geodetic(self.X, self.Y, self.Z)
        self.kn, self.ke, self.kz = cc.vector_cartesian_to_geodetic(self.kvec[:,:,0], self.kvec[:,:,1], self.kvec[:,:,2], self.X, self.Y, self.Z)
        # print self.lat.shape, self.kn.shape
            # self.ke.append(ke)
            # self.kn.append(kn)
            # self.kz.append(kz)



def create_dataset(field, radar):
    # input: Field object, Radar object

    Vx = field.interpVx(np.array([radar.X, radar.Y, radar.Z]).T)
    Vy = field.interpVy(np.array([radar.X, radar.Y, radar.Z]).T)
    Vz = field.interpVz(np.array([radar.X, radar.Y, radar.Z]).T)
    Vvec = np.array([Vx, Vy, Vz]).T

    time0 = (dt.datetime(2019,5,28,0,0)-dt.datetime.utcfromtimestamp(0)).total_seconds()
    times = np.array([[time0+t*60., time0+(t+1)*60.] for t in range(10)])

    Vlos = []
    for t in times:
        Vlos.append(np.einsum('...i,...i->...',radar.kvec, Vvec))
    Vlos = np.array(Vlos)

    print Vvec.shape, radar.kvec.shape, Vlos.shape
    dVlos = np.full(Vlos.shape, 10.)

    # create output hdf5 file
    s = Vlos.shape
    fit_array = np.full((s[0],s[1],s[2],6,4),np.nan)
    fit_array[:,:,:,0,3] = Vlos
    print fit_array.shape

    err_array = np.full((s[0],s[1],s[2],6,4),np.nan)
    err_array[:,:,:,0,3] = dVlos

    ne = np.full(Vlos.shape, 1e11)
    beam_codes = np.array([np.arange(len(radar.elevation)), radar.azimuth, radar.elevation, np.zeros(len(radar.elevation))])


    filename = 'synthetic_data.h5'
    with tables.open_file(filename, mode='w') as file:
        file.create_group('/','FittedParams')
        file.create_group('/','Geomag')
        file.create_group('/','Time')
        file.create_group('/','Site')

        file.create_array('/','BeamCodes',beam_codes)

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