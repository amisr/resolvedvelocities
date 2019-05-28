import numpy as np
from apexpy import Apex
from scipy import interpolate
import coord_convert as cc


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
        self.Vx, self.Vy, self.Vz = cc.vector_geodetic_to_cartesian(self.field[:,0], self.field[:,1], self.field[:,2], self.latitude, self.longitude, self.altitude)

    def create_interpolators(self):

        self.interpVx = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vx)
        self.interpVy = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vy)
        self.interpVz = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vz)


class Radar(object):
    def __init__(self, site, azimuth, elevation, range_step):

        self.X0, self.Y0, self.Z0 = cc.geodetic_to_cartesian(site[0], site[1], site[2])
        self.bin_locations(azimuth, elevation, range_step)

    def bin_locations(self, az, el, rs):
        # a lot of this can probably be handled more elegantly with clever broadcasting, but this works for now

        ranges = np.arange(80.,800., rs)

        self.X = []
        self.Y = []
        self.Z = []
        self.kvec = []

        for a, e in zip(az,el):
            kx = np.cos(e)*np.sin(a)
            ky = np.cos(e)*np.cos(a)
            kz = np.sin(e)

            self.X.append(kx*ranges + self.X0)
            self.Y.append(ky*ranges + self.Y0)
            self.Z.append(kz*ranges + self.Z0)
            self.kvec.append(np.tile(np.array([kx, ky, kz]), (len(ranges),1)))

        self.X = np.array(self.X)
        self.Y = np.array(self.Y)
        self.Z = np.array(self.Z)
        self.kvec = np.array(self.kvec)


