# SyntheticData.py

import numpy as np
import datetime as dt
# from apexpy import Apex
# from scipy import interpolate
import coord_convert as cc
import tables
# import sys
# import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.cm as cm
from matplotlib.colors import Normalize

import configparser

# import sys
# sys.path.append("..")
from ..ResolveVectors import ResolveVectors

class SyntheticData(object):

    def __init__(self, field, radar):
        self.create_dataset(field, radar)
        self.save_dataset(radar)


    def create_dataset(self,field, radar):
        # input: Field object, Radar object

        # interpolate the field to the radar bin locations
        Vx = field.interpVx(np.array([radar.X, radar.Y, radar.Z]).T)
        Vy = field.interpVy(np.array([radar.X, radar.Y, radar.Z]).T)
        Vz = field.interpVz(np.array([radar.X, radar.Y, radar.Z]).T)
        Vvec = np.array([Vx, Vy, Vz]).T

        # create unix time array
        # Time array is 10 time steps (of integration period defined by the radar portion of the config file) after midnight
        #   on January 1 of the year defined in the apex_year portion of the config file.  Because the field is defined manually
        #   and not based on some empirical model, the time really doesn't matter and is mostly included to be consistent
        #   with the read data file format.
        time0 = (dt.datetime(field.apex_year,1,1)-dt.datetime.utcfromtimestamp(0)).total_seconds()
        self.times = np.array([[time0+t*radar.integration_period, time0+(t+1)*radar.integration_period] for t in range(10)])

        # calculate LoS velocity for each bin by taking the dot product of the radar kvector and the interpolated field
        Vlos = np.tile(np.einsum('...i,...i->...',radar.kvec, Vvec), (len(self.times),1,1))
        # assume constant error
        dVlos = np.full(Vlos.shape, radar.vel_error)

        # create fit and error arrays that match the shape of whats in the processed fitted files
        s = Vlos.shape
        self.fit_array = np.full((s[0],s[1],s[2],6,4),np.nan)
        self.fit_array[:,:,:,0,3] = Vlos
        self.err_array = np.full((s[0],s[1],s[2],6,4),np.nan)
        self.err_array[:,:,:,0,3] = dVlos

        self.chi2 = np.full(s, 1.0)
        self.fitcode = np.full(s, 1)

        # generate dummy density array
        self.ne = np.full(s, 1e11)


    def save_dataset(self, radar):

        # create output hdf5 file
        with tables.open_file(radar.output_filename, mode='w') as file:
            file.create_group('/','FittedParams')
            file.create_group('/','Geomag')
            file.create_group('/','Time')
            file.create_group('/','Site')

            file.create_group('/FittedParams', 'FitInfo')

            file.create_array('/','BeamCodes',radar.beam_codes)

            file.create_array('/FittedParams','Fits',self.fit_array)
            file.create_array('/FittedParams','Errors',self.err_array)
            file.create_array('/FittedParams','Ne',self.ne)

            file.create_array('/FittedParams/FitInfo', 'chi2', self.chi2)
            file.create_array('/FittedParams/FitInfo', 'fitcode', self.fitcode)

            file.create_array('/Geomag','Latitude',radar.lat)
            file.create_array('/Geomag','Longitude',radar.lon)
            file.create_array('/Geomag','Altitude',radar.alt*1000.)

            file.create_array('/Geomag','ke',radar.ke)
            file.create_array('/Geomag','kn',radar.kn)
            file.create_array('/Geomag','kz',radar.kz)

            file.create_array('/Time','UnixTime',self.times)

            file.create_array('/Site','Latitude',radar.site_coords[0])
            file.create_array('/Site','Longitude',radar.site_coords[1])
            file.create_array('/Site','Altitude',radar.site_coords[2])


    def eval_vvels(self, config_file):

        rv = ResolveVectors(config_file)
        rv.read_data()
        rv.filter_data()
        rv.transform()
        rv.ion_upflow_correction()
        rv.bin_data()
        rv.get_integration_periods()
        rv.compute_vector_velocity()
        rv.compute_electric_field()
        rv.compute_geodetic_output()

        return rv


    def plot(self, field, radar, rv):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, vx, vy, vz in zip(field.X, field.Y, field.Z, field.Vx, field.Vy, field.Vz):
            ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), linewidths=0.5, color='lightgreen')


        for x, y, z, kx, ky, kz, v in zip(radar.X.ravel(), radar.Y.ravel(), radar.Z.ravel(), radar.kx.ravel(), radar.ky.ravel(), radar.kz.ravel(), self.fit_array[0,:,:,0,3].ravel()):
            ax.quiver(x, y, z, kx*v*1000, ky*v*1000, kz*v*1000, color=quiver_color(v,-500.,500.,'coolwarm'))

        X, Y, Z = cc.geodetic_to_cartesian(rv.bin_glat, rv.bin_glon, rv.bin_galt)
        Vx, Vy, Vz = cc.vector_geodetic_to_cartesian(rv.Velocity_gd[0,:,:,1], rv.Velocity_gd[0,:,:,0], rv.Velocity_gd[0,:,:,2], rv.bin_glat, rv.bin_glon, rv.bin_galt)

        for x, y, z, vx, vy, vz in zip(X.ravel(), Y.ravel(), Z.ravel(), Vx.ravel(), Vy.ravel(), Vz.ravel()):
            ax.quiver(x, y, z, vx, vy, vx, length=0.4*np.sqrt(vx**2+vy**2+vz**2))

        plt.show()

    def compare_components(self, field, rv):

        # convert bin locations to ECEF
        bin_glat, bin_glon, __ = field.A.apex2geo(rv.bin_mlat, rv.bin_mlon, 300.)
        X, Y, Z = cc.geodetic_to_cartesian(bin_glat, bin_glon, np.full(bin_glat.shape,300.))
        # interpolate field to bin locations
        Vx = field.interpVx(np.array([X,Y,Z]).T)
        Vy = field.interpVy(np.array([X,Y,Z]).T)
        Vz = field.interpVz(np.array([X,Y,Z]).T)
        # convert field components to apex
        Vn, Ve, Vu = cc.vector_cartesian_to_geodetic(Vx, Vy, Vz, X, Y, Z)
        f1,f2,f3,g1,g2,g3,d1,d2,d3,e1,e2,e3 = field.A.basevectors_apex(bin_glat, bin_glon, np.full(bin_glat.shape,300.))
        Ve1 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d1.T)
        Ve2 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d2.T)
        Ve3 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d3.T)

        fig = plt.figure(figsize=(10,10))

        ax = fig.add_subplot(311)
        ax.errorbar(np.arange(rv.Velocity.shape[1]),rv.Velocity[0,:,0], yerr=np.sqrt(rv.VelocityCovariance[0,:,0,0]))
        ax.plot(Ve1)
        ax.set_title('Ve1')

        ax = fig.add_subplot(312)
        ax.errorbar(np.arange(rv.Velocity.shape[1]),rv.Velocity[0,:,1], yerr=np.sqrt(rv.VelocityCovariance[0,:,1,1]))
        ax.plot(Ve2)
        ax.set_title('Ve2')

        ax = fig.add_subplot(313)
        ax.errorbar(np.arange(rv.Velocity.shape[1]),rv.Velocity[0,:,2], yerr=np.sqrt(rv.VelocityCovariance[0,:,2,2]))
        ax.plot(Ve3)
        ax.set_title('Ve3')

        plt.show()



def quiver_color(v,vmin,vmax,cmap):
    # get quiver colors (nessisary because 3D quiver routine does wierd stuff with arrow heads)
    norm = Normalize(vmin=vmin,vmax=vmax)
    c = norm(v)
    # c = np.concatenate((c, np.repeat(c, 2)))
    c = getattr(plt.cm,cmap)(c)
    return c    
