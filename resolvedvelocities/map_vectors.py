# map_vectors.py

import ResolveVectors as rv
import numpy as np
from apexpy import Apex
import tables
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import cartopy.crs as ccrs
from matplotlib.colors import Normalize

import coord_convert as cc
import synthetic_data as synth

def plot_raw():
    vvels = rv.ResolveVectors(config='config.ini')

    # get input data
    vvels.read_data()
    vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()

    # idx = 40
    idx = 4

    lat = vvels.lat
    lon = vvels.lon
    alt = vvels.alt/1000.
    ke = vvels.ke
    kn = vvels.kn
    kz = vvels.kz
    if vvels.integration_time:
        # if post integration specified, get vlos for all times in integration period
        vlos = vvels.vlos[vvels.int_idx[idx]]
    else:
        # if no post integration, just get vlos for single index
        vlos = vvels.vlos[idx][None,:]


    x, y, z = cc.geodetic_to_cartesian(lat, lon, alt)


    # calculate vector velocities
    vvels.compute_vectors()
    vvels.compute_geodetic_output()

    lato = vvels.bin_glat
    lono = vvels.bin_glon
    alto = vvels.bin_galt
    vv = vvels.Velocity_gd[idx,:,:]
    vmag = vvels.Vgd_mag[idx]

    fout = np.isfinite(vmag)

    xo, yo, zo = cc.geodetic_to_cartesian(lato, lono, alto)
    vxo, vyo, vzo = cc.vector_geodetic_to_cartesian(vv[:,1],vv[:,0],vv[:,2], lato, lono, alto)


    scale = 100.

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')

    # plot input and output points color coded with bin
    for i,bidx in enumerate(vvels.bin_idx):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.scatter(xo[i], yo[i], zo[i], color=color)
        ax.scatter(x[bidx], y[bidx], z[bidx], color=color)

    # plot input los vectors - loop over differnt post integration records
    for v in vlos:
        fin = np.isfinite(v)
        vx, vy, vz = cc.vector_geodetic_to_cartesian(kn*v, ke*v, kz*v, lat, lon, alt)
        ax.quiver(x[fin], y[fin], z[fin], vx[fin]*scale, vy[fin]*scale, vz[fin]*scale, color=quiver_color(v[fin],-500.,500.,'bwr'))

    # plot output resolved vectors
    ax.quiver(xo[fout], yo[fout], zo[fout], vxo[fout]*scale, vyo[fout]*scale, vzo[fout]*scale, color=quiver_color(vmag[fout],0.,1000.,'Greys'))

    # fig = plt.figure(figsize=(15,10))
    # ax1 = fig.add_subplot(121)
    # ax2 = fig.add_subplot(122,projection=ccrs.LambertConformal())
    # # ax2 = fig.add_subplot(111,projection=ccrs.PlateCarree())
    # ax2.gridlines()

    # ax1.scatter(vvels.lon, vvels.lat)
    # ax1.quiver(vvels.lon, vvels.lat, vvels.ke, vvels.kn, width=0.003)

    # ax2.scatter(vvels.lon, vvels.lat, transform=ccrs.PlateCarree())
    # ax2.quiver(vvels.lon, vvels.lat, vvels.ke, vvels.kn, width=0.003, transform=ccrs.PlateCarree())

    plt.show()

def quiver_color(v,vmin,vmax,cmap):
    # get quiver colors (nessisary because 3D quiver routine does wierd stuff with arrow heads)
    norm = Normalize(vmin=vmin,vmax=vmax)
    c = norm(v)
    c = np.concatenate((c, np.repeat(c, 2)))
    c = getattr(plt.cm,cmap)(c)
    return c    


def plot_mag():

    vvels = rv.ResolveVectors(config='config.ini')

    idx = 4

    # get input data
    vvels.read_data()
    # vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()
    vvels.compute_vectors()

    A = vvels.A
    vlos = vvels.vlos
    vv = vvels.Velocity[idx]

    alt = np.full(vvels.mlat.shape, 300.)
    bin_alt = np.full(vvels.bin_mlat.shape, 300.)


    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')

    s1 = 1.
    s2 = 0.001
    scale = 0.0001

    # plot input and output points color coded with bin
    for mlat,mlon,bidx in zip(vvels.bin_mlat, vvels.bin_mlon, vvels.bin_idx):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.scatter(mlon, mlat, 300., color=color)
        ax.scatter(vvels.mlon[bidx], vvels.mlat[bidx], alt[bidx], color=color)

    # ax.quiver(vvels.mlon, vvels.mlat, alt, A[:,0]*s1, A[:,1]*s1, A[:,2]*s1)
    ax.quiver(vvels.mlon, vvels.mlat, alt, A[:,0]*vlos*s2, A[:,1]*vlos*s2, A[:,2]*vlos*s2)

    ax.quiver(vvels.bin_mlon, vvels.bin_mlat, bin_alt, vv[:,0]*scale, vv[:,1]*scale, vv[:,2]*scale)


    plt.show()

def plot_synth():

    synth_grid = np.meshgrid(np.linspace(62.,70.,10), np.linspace(260.,275.,10))
    velocity = np.tile(np.array([500.,0.,0.]), (10,10,1))
    synth_field = synth.Field(synth_grid[0], synth_grid[1], velocity, 300.)

    az = [14.04,-154.30,-34.69,75.03]
    el = [90.0, 77.5, 66.09, 65.56]
    site = [65.13, -147.47, 0.213]
    radar = synth.Radar(site, az, el, 70.)

    synth.create_dataset(synth_field, radar)


    vvels = rv.ResolveVectors(config='config.ini')

    # get input data
    vvels.read_data()
    # vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()

    idx = 4

    vlos = vvels.vlos[idx]

    x, y, z = cc.geodetic_to_cartesian(vvels.lat, vvels.lon, vvels.alt/1000.)
    kx, ky, kz = cc.vector_geodetic_to_cartesian(vvels.kn, vvels.ke, vvels.kz, vvels.lat, vvels.lon, vvels.alt/1000.)
    vx, vy, vz = cc.vector_geodetic_to_cartesian(vvels.kn*vlos, vvels.ke*vlos, vvels.kz*vlos, vvels.lat, vvels.lon, vvels.alt/1000.)

    # calculate vector velocities
    vvels.compute_vectors()
    vvels.compute_geodetic_output()

    vv = vvels.Velocity_gd[idx,:,:]
    # vmag = vvels.Vgd_mag[idx]

    xo, yo, zo = cc.geodetic_to_cartesian(vvels.bin_glat, vvels.bin_glon, vvels.bin_galt)
    vxo, vyo, vzo = cc.vector_geodetic_to_cartesian(vv[:,1],vv[:,0],vv[:,2], vvels.bin_glat, vvels.bin_glon, vvels.bin_galt)



    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')

    s1=100.
    ax.scatter(synth_field.X,synth_field.Y,synth_field.Z, s=0.1)
    ax.quiver(synth_field.X,synth_field.Y,synth_field.Z, synth_field.Vx*s1, synth_field.Vy*s1, synth_field.Vz*s1)

    s2 = 100000.
    ax.scatter(radar.X, radar.Y, radar.Z)
    # ax.quiver(radar.X, radar.Y, radar.Z, radar.kvec[:,:,0]*s2,radar.kvec[:,:,1]*s2,radar.kvec[:,:,2]*s2, color='orange')

    # ax.scatter(x, y, z)
    # ax.quiver(x, y, z, kx*s2, ky*s2, kz*s2)

    s3 = 1000.
    ax.quiver(x, y, z, vx*s3, vy*s3, vz*s3, color='green')

    ax.scatter(xo, yo, zo)
    ax.quiver(xo, yo, zo, vxo*s1, vyo*s1, vzo*s1)

    plt.show()



def read_vvels_file(filename):
    
    with tables.open_file(filename,mode='r') as file:
        utime = file.get_node('/UnixTime').read()
        utime = (utime[:,0]+utime[:,1])/2.
        # print utime

        # # convert targtime to unix timestamp
        # targtstmp = (time-dt.datetime.utcfromtimestamp(0)).total_seconds()
        # # find index of time in timearray that is closest to targtime
        # t = np.argmin(np.abs(utime-targtstmp))
        # print t
        idx = 10

        E = file.get_node('/ElectricField')[10,:,:]
        VE = file.get_node('/Velocity')[10,:,:]
        # print E.shape, VE.shape

        mlon = file.get_node('/MagneticLongitude').read()
        mlat = file.get_node('/MagneticLatitude').read()
        # print mlat, mlon

    return VE, E, mlat, mlon


def map_vec(alt):
    # calculates plasma drift velocity and electric field vectors in geodetic components at a particular altitude

    VE, E, mlat, mlon = read_vvels_file('test_vvels.h5')
    A = Apex(2019)
    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(mlat,mlon,alt,coords='apex')

    VEgd = (VE[:,0]*e1 + VE[:,1]*e2 + VE[:,2]*e3).T
    Egd = (E[:,0]*d1 + E[:,1]*e2 + E[:,2]*e3).T

    return VEgd, Egd

