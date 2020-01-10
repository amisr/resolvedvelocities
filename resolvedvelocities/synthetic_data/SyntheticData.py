# SyntheticData.py

import numpy as np
import datetime as dt
import pymap3d as pm
import tables
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

import configparser

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
            file.create_array('/Geomag','Altitude',radar.alt)

            file.create_array('/Geomag','ke',radar.ke)
            file.create_array('/Geomag','kn',radar.kn)
            file.create_array('/Geomag','kz',radar.ku)

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
        rv.compute_apex_velocity()
        rv.compute_electric_field()
        rv.compute_geodetic_output()

        return rv


    def plot(self, field, radar, rv):
        # This tends to create plots that show significant disagreement between the "true" field and the resolved field, however,
        #   numerically, both the components of both are very similar.  No sure what causes this difference - possibly some kind
        #   of wierd plotting perspective effect?.

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111,projection='3d')

        for x, y, z, vx, vy, vz in zip(field.X, field.Y, field.Z, field.Vx, field.Vy, field.Vz):
            ax.quiver(x, y, z, vx, vy, vz, length=0.4*np.sqrt(vx**2+vy**2+vz**2), linewidths=0.5, color='lightgreen')


        for x, y, z, kx, ky, kz, v in zip(radar.X.ravel(), radar.Y.ravel(), radar.Z.ravel(), radar.kx.ravel(), radar.ky.ravel(), radar.kz.ravel(), self.fit_array[0,:,:,0,3].ravel()):
            ax.quiver(x, y, z, kx*v*1000, ky*v*1000, kz*v*1000, color=quiver_color(v,-500.,500.,'coolwarm'))

        X, Y, Z = pm.geodetic2ecef(rv.bin_glat, rv.bin_glon, rv.bin_galt*1000.)
        Vx, Vy, Vz = pm.enu2uvw(rv.Velocity_gd[0,:,:,0], rv.Velocity_gd[0,:,:,1], rv.Velocity_gd[0,:,:,2], rv.bin_glat, rv.bin_glon)

        for x, y, z, vx, vy, vz in zip(X.ravel(), Y.ravel(), Z.ravel(), Vx.ravel(), Vy.ravel(), Vz.ravel()):
            ax.quiver(x, y, z, vx, vy, vx, length=0.4*np.sqrt(vx**2+vy**2+vz**2))


        for bn in rv.binvert:
            # interpolate between each bin vertex
            bin_edge = np.empty((0,2))
            for i in range(len(bn)):
                bin_edge = np.concatenate((bin_edge,np.linspace(bn[i-1], bn[i], 10)))
            # convert edge from MARP to cartesian
            glat, glon, _ = rv.marp.marp2geo(bin_edge[:,0], bin_edge[:,1], 300.)
            x, y, z = pm.geodetic2ecef(glat, glon, np.full(glat.shape, 300.)*1000.)
            ax.plot(x, y, z)

        plt.show()

    def check_assumptions(self, field, radar, rv):
        # check the underling assumption that Ve1, Ve2, Ve3 are consistent over the bins
        # WARNING - Mapping is this funciton is deceptive and useless.  This is probably because quivers in the polar cap
        #       aren't being plotted correctly with the chosen projection.

        X = radar.X.flatten()
        Y = radar.Y.flatten()
        Z = radar.Z.flatten()

        glat, glon, galt = pm.ecef2geodetic(X, Y, Z)
        mlat, mlon = rv.marp.geo2marp(glat, glon, galt)

        # field at radar
        Vx = field.interpVx(np.array([X, Y, Z]).T)
        Vy = field.interpVy(np.array([X, Y, Z]).T)
        Vz = field.interpVz(np.array([X, Y, Z]).T)
        Ve, Vn, Vu = pm.uvw2enu(Vx, Vy, Vz, glat, glon)
        Vgd = np.array([Ve, Vn, Vu])

        d1, d2, d3, e1, e2, e3 = rv.marp.basevectors_marp(glat, glon, galt)
        Ve1 = np.einsum('i...,i...->...', d1, Vgd)
        Ve2 = np.einsum('i...,i...->...', d2, Vgd)
        Ve3 = np.einsum('i...,i...->...', d3, Vgd)

        # get base vectors at all points where the field is defined
        d1, d2, d3, e1, e2, e3 = rv.marp.basevectors_marp(field.field_coords[0,:], field.field_coords[1,:],field.field_coords[2,:])

        bins = []
        binsgd = []
        for bn in rv.binvert:
            # interpolate between each bin vertex
            bin_edge = np.empty((0,2))
            for i in range(len(bn)):
                bin_edge = np.concatenate((bin_edge,np.linspace(bn[i-1], bn[i], 10)))
            bins.append(bin_edge)
            # convert edge from MARP to cartesian
            lat, lon, _ = rv.marp.marp2geo(bin_edge[:,0], bin_edge[:,1], 300.)
            binsgd.append(np.array([lat, lon]).T)


        fig = plt.figure(figsize=(17,5))
        gs = gridspec.GridSpec(1,4)
        gs.update(left=0.01,right=0.95)

        # ax = plt.subplot(gs[0,0],projection=ccrs.LambertConformal(central_longitude=np.mean(glon),central_latitude=np.mean(glat)))
        ax = plt.subplot(gs[0,0],projection=ccrs.NorthPolarStereo())
        ax.coastlines()
        ax.scatter(glon, glat, s=2, color='lightgrey', transform=ccrs.PlateCarree())
        ax.quiver(field.field_coords[1,:], field.field_coords[0,:], field.field_values[:,0], field.field_values[:,1], width=0.005, transform=ccrs.PlateCarree())
        # ax.quiver(field.field_coords[1,:], field.field_coords[0,:], e1[0,:], e1[1,:], width=0.005, color='green', label=r'$\hat{e_1}$', transform=ccrs.PlateCarree())
        # ax.quiver(field.field_coords[1,:], field.field_coords[0,:], e2[0,:], e2[1,:], width=0.005, color='red', label=r'$\hat{e_2}$', transform=ccrs.PlateCarree())
        for bin_edge in binsgd:
            ax.plot(bin_edge[:,1],bin_edge[:,0],color='black', transform=ccrs.PlateCarree())
        ax.legend()

        for i, (param,title) in enumerate(zip([Ve1,Ve2,Ve3],['Ve1','Ve2','Ve3'])):
            ax = plt.subplot(gs[0,i+1])
            c = ax.scatter(mlon, mlat, c=param)
            for bin_edge in bins:
                ax.plot(bin_edge[:,1],bin_edge[:,0],color='black')
            ax.set_title(title)
            plt.colorbar(c)

        plt.show()
        # plt.savefig('synth_data_assumptions.png')

    def check_components(self, field, rv):

        targalt = 300.
        aidx = np.argmin(np.abs(rv.bin_galt[:,0]-targalt))

        # convert bin locations to ECEF
        X, Y, Z = pm.geodetic2ecef(rv.bin_glat[aidx,:], rv.bin_glon[aidx,:], rv.bin_galt[aidx,:]*1000.)
        # interpolate field to bin locations
        Vx = field.interpVx(np.array([X,Y,Z]).T)
        Vy = field.interpVy(np.array([X,Y,Z]).T)
        Vz = field.interpVz(np.array([X,Y,Z]).T)
        # convert field components to apex
        Ve, Vn, Vu = pm.uvw2enu(Vx, Vy, Vz, rv.bin_glat[aidx,:], rv.bin_glon[aidx,:])
        d1,d2,d3,e1,e2,e3 = rv.marp.basevectors_marp(rv.bin_glat[aidx,:], rv.bin_glon[aidx,:], rv.bin_galt[aidx,:])
        Ve1 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d1.T)
        Ve2 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d2.T)
        Ve3 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d3.T)

        # get bin outlines
        bins = []
        binsgd = []
        for bn in rv.binvert:
            # interpolate between each bin vertex
            bin_edge = np.empty((0,2))
            for i in range(len(bn)):
                bin_edge = np.concatenate((bin_edge,np.linspace(bn[i-1], bn[i], 10)))
            bins.append(bin_edge)
            # convert edge from MARP to cartesian
            lat, lon, _ = rv.marp.marp2geo(bin_edge[:,0], bin_edge[:,1], 300.)
            binsgd.append(np.array([lat, lon]).T)


        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(4,2)
        gs.update(bottom=0.05, top=0.95, left=0.1, right=0.95, hspace=0.3)

        ax = plt.subplot(gs[0,0])
        ax.quiver(rv.bin_mlon,rv.bin_mlat,Ve1,Ve2, scale=10000, color='blue')
        ax.quiver(rv.bin_mlon,rv.bin_mlat,rv.Velocity[0,:,0],rv.Velocity[0,:,1], scale=10000, color='orange')
        for bin_edge in bins:
            ax.plot(bin_edge[:,1],bin_edge[:,0],color='black')
        # ax.set_xlim([1.5*min(rv.bin_mlon)-0.5*max(rv.bin_mlon),1.5*max(rv.bin_mlon)-0.5*min(rv.bin_mlon)])
        # ax.set_ylim([1.5*min(rv.bin_mlat)-0.5*max(rv.bin_mlat),1.5*max(rv.bin_mlat)-0.5*min(rv.bin_mlat)])
        ax.set_title('Magnetic/Native')

        for i, (param,title) in enumerate(zip([Ve1,Ve2,Ve3],['Ve1','Ve2','Ve3'])):
            ax = plt.subplot(gs[i+1,0])
            ax.plot(param, color='blue')
            ax.errorbar(np.arange(rv.Velocity.shape[1]),rv.Velocity[0,:,i], yerr=np.sqrt(rv.VelocityCovariance[0,:,i,i]), color='orange')
            ax.set_title(title)

        ax = plt.subplot(gs[0,1],projection=ccrs.LambertConformal(central_longitude=np.mean(field.field_coords[1,:]),central_latitude=np.mean(field.field_coords[0,:])))
        ax.coastlines()

        ax.quiver(rv.bin_glon[aidx,:], rv.bin_glat[aidx,:], Ve, Vn, scale=2000, color='blue', transform=ccrs.PlateCarree())
        ax.quiver(rv.bin_glon[aidx,:], rv.bin_glat[aidx,:], rv.Velocity_gd[0,aidx,:,0], rv.Velocity_gd[0,aidx,:,1], scale=2000, color='orange', transform=ccrs.PlateCarree())
        for bin_edge in binsgd:
            ax.plot(bin_edge[:,1],bin_edge[:,0],color='black', transform=ccrs.PlateCarree())
        # ax.set_extent([min(rv.bin_glon[aidx,:]),max(rv.bin_glon[aidx,:]),min(rv.bin_glat[aidx,:]),max(rv.bin_glat[aidx,:])])
        # ax.set_extent([1.5*min(rv.bin_glon[aidx,:])-0.5*max(rv.bin_glon[aidx,:]),1.5*max(rv.bin_glon[aidx,:])-0.5*min(rv.bin_glon[aidx,:]),1.5*min(rv.bin_glat[aidx,:])-0.5*max(rv.bin_glat[aidx,:]),1.5*max(rv.bin_glat[aidx,:])-0.5*min(rv.bin_glat[aidx,:])])
        ax.set_title('Geodetic')

        for i, (param,title) in enumerate(zip([Ve,Vn,Vu],['Ve','Vn','Vu'])):
            ax = plt.subplot(gs[i+1,1])
            ax.plot(param, color='blue')
            ax.errorbar(np.arange(rv.Velocity_gd.shape[2]),rv.Velocity_gd[0,aidx,:,i], yerr=np.sqrt(rv.VelocityCovariance_gd[0,aidx,:,i,i]),color='orange')
            ax.set_title(title)

        plt.show()
        # plt.savefig('synth_data_components.png')



def quiver_color(v,vmin,vmax,cmap):
    # get quiver colors (nessisary because 3D quiver routine does wierd stuff with arrow heads)
    norm = Normalize(vmin=vmin,vmax=vmax)
    c = norm(v)
    # c = np.concatenate((c, np.repeat(c, 2)))
    c = getattr(plt.cm,cmap)(c)
    return c
