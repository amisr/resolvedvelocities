# Radar.py
import numpy as np
# import datetime as dt
# from apexpy import Apex
# from scipy import interpolate
import coord_convert as cc
# import tables
# import sys
# import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.cm as cm
# from matplotlib.colors import Normalize

import configparser

class Radar(object):

    def __init__(self, config):

        self.read_config(config)

        bc = np.loadtxt(self.beamcode_filename)
        idx = np.where(np.in1d(bc[:,0],self.beamcodes))[0]
        self.beam_codes = bc[idx,:]

        self.X0, self.Y0, self.Z0 = cc.geodetic_to_cartesian(self.site_coords[0], self.site_coords[1], self.site_coords[2])
        self.get_gate_locations(self.beam_codes[:,1], self.beam_codes[:,2], float(self.range_step))
        self.geodetic_locations()
        # self.plot_radar()

    def read_config(self, config_file):
        config = configparser.ConfigParser()
        config.read(config_file)
        self.__dict__.update(config.items('RADAR'))
        self.__dict__.update((name,eval(value)) for name, value in self.__dict__.items())



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
        kx, ky, kz = cc.vector_geodetic_to_cartesian(kn, ke, ku, self.site_coords[0], self.site_coords[1], self.site_coords[2])

        # calculate position of each range gate in ECEF
        self.X = kx[:,None]*ranges + self.X0
        self.Y = ky[:,None]*ranges + self.Y0
        self.Z = kz[:,None]*ranges + self.Z0

        # form array of k vector at each range gate
        self.kvec = np.array([np.tile(np.array([x,y,z]), (len(ranges),1)) for x, y, z in zip(kx,ky,kz)])

        self.kx = np.tile(kx, (len(ranges),1)).T
        self.ky = np.tile(ky, (len(ranges),1)).T
        self.kz = np.tile(kz, (len(ranges),1)).T


    def geodetic_locations(self):
        # calculate gate position and k vectors in geodetic coordinates
        self.lat, self.lon, self.alt = cc.cartesian_to_geodetic(self.X, self.Y, self.Z)
        self.kn, self.ke, self.kz = cc.vector_cartesian_to_geodetic(self.kx, self.ky, self.kz, self.X, self.Y, self.Z)

    def plot_radar(self):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, kx, ky, kz in zip(self.X, self.Y, self.Z, self.kx, self.ky, self.kz):
            ax.quiver(x, y, z, kx, ky, kz, length=100000, color='orange')


        plt.show()
