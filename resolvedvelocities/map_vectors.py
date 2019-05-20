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

def plot_raw():
    vvels = rv.ResolveVectors()

    # get input data
    vvels.read_data()
    vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()

    idx = 40

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

    lato = vvels.gdlat
    lono = vvels.gdlon
    alto = vvels.gdalt
    vv = vvels.Velocity_gd[idx,:,:]
    vmag = vvels.Vgd_mag[idx]

    fout = np.isfinite(vmag)

    xo, yo, zo = cc.geodetic_to_cartesian(lato, lono, alto)
    vxo, vyo, vzo = cc.vector_geodetic_to_cartesian(vv[:,0],vv[:,1],vv[:,2], lato, lono, alto)


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

    vvels = rv.ResolveVectors()

    idx = 40

    # get input data
    vvels.read_data()
    vvels.filter_data()
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()
    vvels.compute_vectors()

    A = vvels.A
    vv = vvels.Velocity[idx]


    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)

    scale = 0.001

    # plot input and output points color coded with bin
    for mlat,mlon,bidx in zip(vvels.bin_mlat, vvels.bin_mlon, vvels.bin_idx):
        color = next(ax._get_lines.prop_cycler)['color']
        ax.scatter(mlon, mlat, color=color)
        ax.scatter(vvels.mlon[bidx], vvels.mlat[bidx], color=color)

    ax.quiver(vvels.mlon, vvels.mlat, A[:,0]*scale, A[:,1]*scale)

    ax.quiver(vvels.bin_mlon, vvels.bin_mlat, vv[:,0]*scale, vv[:,1]*scale)


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

