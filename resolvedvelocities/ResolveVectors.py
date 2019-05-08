# ResolveVectors.py


class ResolveVectors(object):
    def __init__(self):
        # read config file
        # set class attributes according to config file
        pass


    def read_data(self):
        # read data from standard AMISR fit files
        pass

    def transform(self):
        # transform k vectors from geodetic to geomagnetic
        # use apexpy
        pass

    def bin_data(self):
        # divide data into an arbitrary number of bins
        # bins defined in some way by initial config file
        # each bin has a specified MLAT/MLON
        pass

    def compute_vectors(self):
        # use Heinselman and Nicolls Bayesian reconstruction algorithm to get full vectors
        pass

    def compute_geodetic_output(self):
        # map velocity and electric field to get an array at different altitudes
        pass

    def save_output(self):
        # save output file
        pass