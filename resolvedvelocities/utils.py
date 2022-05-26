# utils.py
# Convenient utility functions needed for resolving vector velocities.

import numpy as np
import datetime as dt
from apexpy import Apex
import tables
import os
import importlib
import socket
import getpass
import platform



def vvels(vlos, dvlos, A, cov, minnumpoints=1):
    # Bayesian inference method described in Heinselman and Nicolls 2008
    # vector velocity algorithm

    # Get indices for finite valued data points
    finite = np.isfinite(vlos)
    num_points = np.sum(finite)
    dof = num_points - 3 # solving for 3 components

    # Filter inputs to only use finite valued data
    vlos = vlos[finite]
    dvlos = dvlos[finite]
    A = A[finite]

    SigmaE = np.diagflat(dvlos**2)
    SigmaV = np.diagflat(cov)

    try:
        # measurement errors and a priori covariance, terms in the inverse
        # (Heinselman and Nicolls 2008 eqn 12)
        # I = (A*SigV*A.T + SigE)^-1
        I = np.linalg.inv(np.einsum('jk,kl,ml->jm',A,SigmaV,A) + SigmaE)
        # calculate velocity estimate (Heinselman and Nicolls 2008 eqn 12)
        V = np.einsum('jk,lk,lm,m->j',SigmaV,A,I,vlos)
        # calculate covariance of velocity estimate
        # (Heinselman and Nicolls 2008 eqn 13)
        term1 = np.einsum('kj,kl,lm->jm',A,np.linalg.inv(SigmaE),A)
        term2 = np.linalg.inv(SigmaV)
        SigV = np.linalg.inv(term1 + term2)

        # chi-squared is meaningless for an underdetermined problem
        if dof < 1:
            chi2 = np.nan
        else:
            model = np.einsum('...i,i->...',A,V)
            chi2 = np.sum((vlos - model)**2 / dvlos**2) / dof

    except np.linalg.LinAlgError:
        V = np.full(3,np.nan)
        SigV = np.full((3,3),np.nan)
        chi2 = np.nan
        num_points = 0

    return V, SigV, chi2, num_points


def parse_list(s):
    return [float(i) for i in s.split(',')]


def lin_interp(x, xp, fp, dfp):
    # Piecewise linear interpolation routine that returns interpolated values
    # and their errors

    # find the indicies of xp that bound each value in x
    # Note: where x is out of range of xp, -1 is used as a place holder
    #   This provides a valid "dummy" index for the array calculations and can
    # be used to identify values to nan in final output
    xpmin = np.nanmin(xp)
    xpmax = np.nanmax(xp)
    i = np.array([np.argwhere((xi>=xp[:-1]) & (xi<xp[1:])).flatten()[0]
                  if ((xi>=xpmin) & (xi<xpmax)) else -1 for xi in x])
    # calculate X
    X = (x-xp[i])/(xp[i+1]-xp[i])
    # calculate interpolated values
    f = (1-X)*fp[i] + X*fp[i+1]
    # calculate interpolation error
    df = np.sqrt((1-X)**2*dfp[i]**2 + X**2*dfp[i+1]**2)
    # replace out-of-range values with NaN
    f[i<0] = np.nan
    df[i<0] = np.nan

    return f, df

def magnitude_direction(A,Sig,e):
    # Calculate the magnitude of vector A and the clockwise angle between
    # vectors e and A
    # Also calculates corresponding errors
    # A = vector
    # Sig = covariance matrix for A
    # e = vector to take the direction relative to
    # ep = e x z (vector perpendicular to e and up)
    # This is all done with an error analysis using addition in quadrature.
    # All calculations are done with matrix algebra using einsum to prevent
    # nested for loops.
    # Input vectors are assumed to have orthogonal components

    # Calculate the Magnitude of input A

    # Some helper matricies
    # dot product of A and A
    AA = np.einsum('...i,...i->...', A, A)
    # matrix multipy A*Sig*transpose(A)
    ASA = np.einsum('...i,...ij,...j->...', A, Sig, A)

    # calculate magnitude and magnitude error
    magnitude = np.sqrt(AA)
    mag_err = np.sqrt(ASA/AA)


    # Now find the angle clockwise around geodetic up in the horizontal plane
    # that is perpendicular to geodetic up, where angle=0 is along the
    # projection 'e' in the horizontal plane

    # Some helper matricies
    # dot product of e and e
    ee = np.einsum('...i,...i->...', e, e)
    # dot product of e and A
    eA = np.einsum('...i,...i->...', e, A)
    # find ep, perpendicular to both e and geodetic up
    ep = np.cross(e,np.array([0,0,1]))
    epep = np.einsum('...i,...i->...', ep, ep)
    epA = np.einsum('...i,...i->...', ep, A)

    # B = ep(e*A)-e(ep*A) = A x (ep x e)
    B = np.einsum('...ij,...i->...ij',ep,eA)-np.einsum('...ij,...i->...ij',e,epA)
    # matrix multipy B*Sig*B (covariance propagation)
    BSB = np.einsum('...i,...ij,...j->...', B, Sig, B)

    # calculate direction and direction error
    direction = np.arctan2(np.sqrt(ee)*epA,np.sqrt(epep)*eA)
    dir_err = np.sqrt(epep*ee*BSB)/(ee*epA**2-epep*eA**2)

    return magnitude, mag_err, direction*180./np.pi, dir_err*180./np.pi


def get_integration_periods(utime, integration_time):
    # calculate new integration periods

    integration_periods = list()
    start_time = None
    num_times = len(utime)
    for i,time_pair in enumerate(utime):

        temp_start_time, temp_end_time = time_pair
        if start_time is None:
            start_time = temp_start_time
        time_diff = temp_end_time - start_time

        if (time_diff >= integration_time) or (i == num_times -1):
            integration_periods.append([start_time, temp_end_time])
            start_time = None
            continue

    return np.array(integration_periods)


def create_time_arrays(integration_periods, site):
    # generate time arrays in various formats
    time_array = np.array([[dt.datetime.utcfromtimestamp(t[0]), dt.datetime.utcfromtimestamp(t[1])] for t in integration_periods])
    year = np.array([[t[0].year, t[1].year] for t in time_array])
    month = np.array([[t[0].month, t[1].month] for t in time_array])
    day = np.array([[t[0].day, t[1].day] for t in time_array])
    doy = np.array([[t[0].timetuple().tm_yday, t[1].timetuple().tm_yday] for t in time_array])
    dtime = np.array([[(t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.), (t[0]-t[0].replace(hour=0,minute=0,second=0)).total_seconds()/(60.*60.)] for t in time_array])

    A = Apex(date=time_array[0,0])
    mlat, mlon = A.geo2apex(site[0], site[1], site[2])
    mlt = np.array([[A.mlon2mlt(mlon,t[0]), A.mlon2mlt(mlon,t[1])] for t in time_array])

    return year, month, day, doy, dtime, mlt


def save_carray(h5, node, data, attributes):
    # Save an array to a hdf5 file as a CArray, which allows built-in compression
    group, name = os.path.split(node)
    atom = tables.Atom.from_dtype(data.dtype)
    arr = h5.create_carray(group, name, atom, data.shape)
    arr[:] = data
    for k, v in attributes.items():
        h5.set_node_attr(node, k, v)

def get_example_config():

    example_config = importlib.resources.read_text(__package__, 'example_config.ini')

    return example_config


def get_processing_info(configfile):

    # Save computer information and config file
    # Computer information:
    computer_info = {}
    computer_info['PythonVersion'] = platform.python_version()
    computer_info['Type'] = platform.machine()
    computer_info['System'] = '{} {} {}'.format(platform.system(),platform.release(),platform.version())
    computer_info['User'] = getpass.getuser()
    hostname = platform.node()
    if len(hostname) == 0:
        hostname = socket.gethostname()
    computer_info['Host'] = hostname

    config_info = {}
    config_info['Path'] = os.path.dirname(os.path.abspath(configfile))
    config_info['Name'] = os.path.basename(configfile)
    with open(configfile, 'r') as f:
        config_info['Contents'] = ''.join(f.readlines())

    software_version = importlib.metadata.version(__package__)

    return computer_info, config_info, software_version
