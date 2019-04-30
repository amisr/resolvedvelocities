# map_vectors.py

from apexpy import Apex
import tables

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

