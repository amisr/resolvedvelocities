# plot.py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import datetime as dt

def plot_components(utime, mlat, vector, covariance):

    time = np.array([dt.datetime.utcfromtimestamp(np.mean(t)) for t in utime])

    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(2,3)

    # TODO: These kinds of plots don't really make much sense if we're allowing arbitrary bin locations
    for i in range(3):
        ax = plt.subplot(gs[0,i])
        ax.pcolormesh(time, mlat, vector[:,:,i].T, vmin=-1500, vmax=1500)

    # for i in range(3):
        ax = plt.subplot(gs[1,i])
        ax.pcolormesh(time, mlat, np.sqrt(covariance[:,:,i,i]).T, vmin=0, vmax=350)

    plt.show()