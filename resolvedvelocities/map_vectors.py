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

    idx = 20

    # form arrays of only finite values (quiver plotting doesn't handle NANs well)
    lat = vvels.lat[np.isfinite(vvels.vlos[idx])]
    lon = vvels.lon[np.isfinite(vvels.vlos[idx])]
    alt = vvels.alt[np.isfinite(vvels.vlos[idx])]/1000.
    ke = vvels.ke[np.isfinite(vvels.vlos[idx])]
    kn = vvels.kn[np.isfinite(vvels.vlos[idx])]
    kz = vvels.kz[np.isfinite(vvels.vlos[idx])]
    vlos = vvels.vlos[idx][np.isfinite(vvels.vlos[idx])]

    x, y, z = cc.geodetic_to_cartesian(lat, lon, alt)
    vx, vy, vz = cc.vector_geodetic_to_cartesian(kn*vlos, ke*vlos, kz*vlos, lat, lon, alt)

    # calculate vector velocities
    vvels.transform()
    vvels.bin_data()
    vvels.get_integration_periods()
    vvels.compute_vectors()
    vvels.compute_geodetic_output()

    # form arrays of only finite values
    lato = vvels.gdlat[np.isfinite(vvels.Velocity_gd[idx,:,0])]
    lono = vvels.gdlon[np.isfinite(vvels.Velocity_gd[idx,:,0])]
    alto = vvels.gdalt[np.isfinite(vvels.Velocity_gd[idx,:,0])]
    vv = vvels.Velocity_gd[idx,:,:]
    vv = vv[np.isfinite(vv[:,0]),:]

    xo, yo, zo = cc.geodetic_to_cartesian(lato, lono, alto)
    vxo, vyo, vzo = cc.vector_geodetic_to_cartesian(vv[:,0],vv[:,1],vv[:,2], lato, lono, alto)



    # get quiver colors (nessisary because 3D quiver routine does wierd stuff with arrow heads)
    norm = Normalize(vmin=-500.,vmax=500.)
    c = norm(vlos)
    c = np.concatenate((c, np.repeat(c, 2)))
    c = plt.cm.bwr(c)

    scale = 100.

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(x, y, z)
    ax.quiver(x, y, z, vx*scale, vy*scale, vz*scale, color=c)

    ax.scatter(xo, yo, zo)
    ax.quiver(xo, yo, zo, vxo*scale, vyo*scale, vzo*scale)
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

# check ke, kn, kz
# in ECEF, bin location minus radar location
# kvec for every bin in beam should be the same

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

