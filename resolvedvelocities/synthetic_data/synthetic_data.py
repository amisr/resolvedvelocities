import numpy as np
import datetime as dt
from apexpy import Apex
from scipy import interpolate
import coord_convert as cc
import tables
import sys
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize

import sys
sys.path.append("..")
from ResolveVectors import ResolveVectors


class Field(object):
    def __init__(self, lat, lon, field, alt):

        self.map_velocity_field(lat.flatten(), lon.flatten(), field.reshape(-1,field.shape[-1]), alt)
        self.convert_to_ECEF()
        self.create_interpolators()
        self.plot_ionosphere()

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


    def plot_ionosphere(self):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, vx, vy, vz in zip(self.X, self.Y, self.Z, self.Vx, self.Vy, self.Vz):
            ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), color='green')

        plt.show()



class Radar(object):
    def __init__(self, site, beams=None, azimuth=None, elevation=None, range_step=50.):
        # beam - list of beam codes
        # azimuth - list of azimuth angle for each beam (in degrees)
        # elevation - list of elevation angle for each beam (in degrees)
        # range_step - step between each range gate (in km)

        amisr_sites = {'PFISR':[65.13,-147.47,0.213], 'RISRN':[74.72955,-94.90576,0.145], 'RISRC':[74.72955,-94.90576,0.145]}

        beam_code_file = 'bcotable_{}.txt'.format(site.lower())

        self.site = amisr_sites[site]
        try:
            # bc = np.loadtxt(os.path.join(os.path.dirname(__file__), 'bcotable.txt'))
            bc = np.loadtxt(beam_code_file)
            idx = np.where(np.in1d(bc[:,0],beams))[0]
            self.beam_codes = bc[idx,:]
        except:
            self.beam_codes = np.array([range(len(azimuth)),azimuth,elevation,np.full(len(azimuth),np.nan)]).T

        self.range_step = range_step
        self.X0, self.Y0, self.Z0 = cc.geodetic_to_cartesian(self.site[0], self.site[1], self.site[2])
        self.get_gate_locations(self.beam_codes[:,1], self.beam_codes[:,2], range_step)
        self.geodetic_locations()

        self.plot_radar()

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

        # print(self.kx.shape, self.X.shape)

        # form array of k vector at each range gate
        self.kvec = np.array([np.tile(np.array([x,y,z]), (len(ranges),1)) for x, y, z in zip(kx,ky,kz)])

        self.kx = np.tile(kx, (len(ranges),1)).T
        self.ky = np.tile(ky, (len(ranges),1)).T
        self.kz = np.tile(kz, (len(ranges),1)).T
        # print(self.kx.shape, self.X.shape)

    def geodetic_locations(self):
        # calculate gate position and k vectors in geodetic coordinates
        self.lat, self.lon, self.alt = cc.cartesian_to_geodetic(self.X, self.Y, self.Z)
        # self.kn, self.ke, self.kz = cc.vector_cartesian_to_geodetic(self.kvec[:,:,0], self.kvec[:,:,1], self.kvec[:,:,2], self.X, self.Y, self.Z)
        self.kn, self.ke, self.kz = cc.vector_cartesian_to_geodetic(self.kx, self.ky, self.kz, self.X, self.Y, self.Z)

    def plot_radar(self):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, kx, ky, kz in zip(self.X, self.Y, self.Z, self.kx, self.ky, self.kz):
            ax.quiver(x, y, z, kx, ky, kz, length=100000, color='orange')


        plt.show()



class SyntheticData(object):

    def __init__(self, field, radar):
        self.create_dataset(field, radar)
        self.eval_vvels()
        self.plot(field, radar)

    def create_dataset(self,field, radar):
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

        self.Vlos = Vlos

        # create fit and error arrays that match the shape of whats in the processed fitted files
        s = Vlos.shape
        fit_array = np.full((s[0],s[1],s[2],6,4),np.nan)
        fit_array[:,:,:,0,3] = Vlos
        err_array = np.full((s[0],s[1],s[2],6,4),np.nan)
        err_array[:,:,:,0,3] = dVlos

        chi2 = np.full(s, 1.0)
        fitcode = np.full(s, 1)

        # generate dummy density array
        ne = np.full(Vlos.shape, 1e11)

        # create output hdf5 file
        filename = 'synthetic_data.h5'
        with tables.open_file(filename, mode='w') as file:
            file.create_group('/','FittedParams')
            file.create_group('/','Geomag')
            file.create_group('/','Time')
            file.create_group('/','Site')

            file.create_group('/FittedParams', 'FitInfo')

            file.create_array('/','BeamCodes',radar.beam_codes)

            file.create_array('/FittedParams','Fits',fit_array)
            file.create_array('/FittedParams','Errors',err_array)
            file.create_array('/FittedParams','Ne',ne)

            file.create_array('/FittedParams/FitInfo', 'chi2', chi2)
            file.create_array('/FittedParams/FitInfo', 'fitcode', fitcode)

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


    def eval_vvels(self):

        rv = ResolveVectors('config.ini')
        rv.read_data()
        rv.filter_data()
        rv.transform()
        rv.ion_upflow_correction()
        rv.bin_data()
        rv.get_integration_periods()
        rv.compute_vector_velocity()
        rv.compute_electric_field()
        rv.compute_geodetic_output()

        print(rv.Velocity_gd.shape)
        print(rv.bin_glat.shape, rv.bin_glon.shape, rv.bin_galt.shape)

        self.X, self.Y, self.Z = cc.geodetic_to_cartesian(rv.bin_glat, rv.bin_glon, rv.bin_galt)
        self.Vx, self.Vy, self.Vz = cc.vector_geodetic_to_cartesian(rv.Velocity_gd[0,:,:,1], rv.Velocity_gd[0,:,:,0], rv.Velocity_gd[0,:,:,2], rv.bin_galt, rv.bin_glon, rv.bin_galt)

    def plot(self, field, radar):

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, vx, vy, vz in zip(field.X, field.Y, field.Z, field.Vx, field.Vy, field.Vz):
            ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), linewidths=0.5, color='lightgreen')


        for x, y, z, kx, ky, kz, v in zip(radar.X.ravel(), radar.Y.ravel(), radar.Z.ravel(), radar.kx.ravel(), radar.ky.ravel(), radar.kz.ravel(), self.Vlos[0].ravel()):
            ax.quiver(x, y, z, kx*v*1000, ky*v*1000, kz*v*1000, color=quiver_color(v,-500.,500.,'coolwarm'))

        for x, y, z, vx, vy, vz in zip(self.X.ravel(), self.Y.ravel(), self.Z.ravel(), self.Vx.ravel(), self.Vy.ravel(), self.Vz.ravel()):
            ax.quiver(x, y, z, vx, vy, vx, length=0.4*np.sqrt(vx**2+vy**2+vz**2))

        plt.show()



def quiver_color(v,vmin,vmax,cmap):
    # get quiver colors (nessisary because 3D quiver routine does wierd stuff with arrow heads)
    norm = Normalize(vmin=vmin,vmax=vmax)
    c = norm(v)
    # c = np.concatenate((c, np.repeat(c, 2)))
    c = getattr(plt.cm,cmap)(c)
    return c    
