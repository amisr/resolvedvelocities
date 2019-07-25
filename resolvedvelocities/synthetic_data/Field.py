# Field.py
import numpy as np
from apexpy import Apex
from scipy import interpolate
import coord_convert as cc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import configparser

class Field(object):

    def __init__(self, config):

        self.read_config(config)

        # convert list from config file to arrays
        self.field_coords = np.array(self.field_coords)
        self.field_values = np.array(self.field_values)

        # initialize Apex object
        self.A = Apex(date=self.apex_year)

        self.map_velocity_field(self.field_coords, self.field_values)
        self.convert_to_ECEF()
        self.create_interpolators()
        # self.plot_ionosphere()

    def read_config(self, config_file):

        config = configparser.ConfigParser()
        config.read(config_file)
        self.__dict__.update(config.items('FIELD'))
        self.__dict__.update((name,eval(value)) for name, value in self.__dict__.items())


    def map_velocity_field(self, coords, field):
        # coords - array (N,3) of geodetic lat, lon, alt
        # field - array (N,3) of geodetic E, N, U components of the velocity field at each position

        # define output altitudes
        altitude = np.arange(50., 1000., 50.)
        # create array in proper shape to be applied to every input coordinate
        self.altitude = np.repeat(altitude,coords.shape[-1])

        # map to diffent altitudes manually - the current expected input/output arrays of apexpy.map_to_height makes this function difficult to use for this purpose
        alat, alon = self.A.geo2apex(coords[0], coords[1], coords[2])
        # find positions at each altitude
        self.latitude, self.longitude, __ = self.A.apex2geo(np.tile(alat,len(altitude)), np.tile(alon,len(altitude)), self.altitude)

        # map field to each altitude
        f = np.array([self.A.map_V_to_height(alat, alon, coords[2], a, field.T).T for a in altitude])
        self.field = f.reshape(-1,f.shape[-1])


    def convert_to_ECEF(self):

        self.X, self.Y, self.Z = cc.geodetic_to_cartesian(self.latitude, self.longitude, self.altitude)
        self.Vx, self.Vy, self.Vz = cc.vector_geodetic_to_cartesian(self.field[:,1], self.field[:,0], self.field[:,2], self.latitude, self.longitude, self.altitude)

    def create_interpolators(self):

        self.interpVx = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vx)
        self.interpVy = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vy)
        self.interpVz = interpolate.LinearNDInterpolator(np.array([self.X, self.Y, self.Z]).T, self.Vz)


    def plot_ionosphere(self):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, vx, vy, vz in zip(self.X, self.Y, self.Z, self.Vx, self.Vy, self.Vz):
            ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), color='green')

        plt.show()
