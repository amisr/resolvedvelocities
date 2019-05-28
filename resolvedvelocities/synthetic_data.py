import numpy as np
from apexpy import Apex


def map_velocity_field(grid, field, alt_in):
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
        field_out.append(full_field.T.reshape((fs[0],fs[1],fs[2])))

        lat, lon, __ = A.apex2geo(grid[0],grid[1],a)
        lat_out.append(lat)
        lon_out.append(lon)
        alt_out.append(np.full(lat.shape,a))

    field_out = np.array(field_out)
    lat_out = np.array(lat_out)
    lon_out = np.array(lon_out)
    alt_out = np.array(alt_out)

    return lat_out, lon_out, alt_out, field_out

