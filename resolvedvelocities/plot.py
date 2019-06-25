# plot.py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import datetime as dt

def plot_components(utime, mlat, mlon, vector, covariance, param='V', titles=None, clim=None, cmap=None):

    # get x-axis (time) tick locations and labels
    time, xticks = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    binloc = ['({} N,\n {} E)'.format(lat, lon) for lat, lon in zip(mlat, mlon)]

    # set defaults
    defaults = {'V': {'titles':[['Ve1 (m/s)', 'Ve2 (m/s)', 'Ve3 (m/s) x 10'],
                                ['errVe1 (m/s)', 'errVe2 (m/s)', 'errVe3 (m/s) x 10']],
                      'clim':[[-1500., 1500.], [0., 350.]],
                      'cmap':['coolwarm', 'viridis']},
                'E': {'titles':[['Ed1 (mV/m)', 'Ed2 (mV/m)', 'Ed3 (mV/m) x 10'],
                                ['errEd1 (mV/m)', 'errEd2 (mV/m)', 'errEd3 (mV/m) x 10']],
                      'clim':[[-75., 75.], [0., 15.]],
                      'cmap':['coolwarm', 'viridis']}}

    # if titles, color limits, and color maps are not specified in function call, use defaults
    if not titles:
        titles = defaults[param]['titles']
    if not clim:
        clim = defaults[param]['clim']
    if not cmap:
        cmap = defaults[param]['cmap']

    # for electric field, multiply values by 1000 to get units of mV/m instead of V/m
    if param=='E':
        vector = vector*1000.
        covariance = covariance*1000.*1000.

    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2,3)
    gs.update(wspace=0.3,hspace=0.2,left=0.05,right=0.9,bottom=0.05,top=0.95)

    for i, A in enumerate([vector,np.sqrt(np.diagonal(covariance,axis1=-1,axis2=-2))]):
        for j in range(3):
            # for along B component, multiply by 10
            if j==2:
                fac=10.
            else:
                fac=1.

            ax = plt.subplot(gs[i,j])
            f = ax.pcolormesh(A[:,:,j].T*fac, vmin=clim[i][0], vmax=clim[i][1], cmap=plt.get_cmap(cmap[i]))
            ax.set_xticks(xticks)
            ax.set_xticklabels(time)
            ax.set_yticks(yticks)
            ax.set_yticklabels(binloc)
            ax.set_title(titles[i][j])
            ax.tick_params(labelsize=8)

        pos = ax.get_position()
        cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
        cbar = fig.colorbar(f, cax=cax)

    plt.savefig('{}components.png'.format(param))
    plt.show()


def plot_magnitude(utime, mlat, mlon, vmag, dvmag, vdir, dvdir, param='V', titles=None, clim=None, cmap=None):

    # get x-axis (time) tick locations and labels
    time, xticks = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    binloc = ['({} N, {} E)'.format(lat, lon) for lat, lon in zip(mlat, mlon)]

    # set defaults
    defaults = {'V': {'titles':['V magnitude (m/s)', 'V magnitude error (m/s)', 'V direction (deg)', 'V direction error (deg)', ''],
                      'clim':[[0.,1500.],[0., 350.],[-180., 180.],[0., 35.]],
                      'cmap':['viridis', 'viridis', 'hsv', 'viridis']},
                'E': {'titles':['E magnitude (mV/m)', 'V magnitude error (mV/m)', 'E direction (deg)', 'E direction error (deg)', ''],
                      'clim':[[0.,75.],[0., 15.],[-180., 180.],[0., 35.]],
                      'cmap':['viridis', 'viridis', 'hsv', 'viridis']}}

    # if titles, color limits, and color maps are not specified in function call, use defaults
    if not titles:
        titles = defaults[param]['titles']
    if not clim:
        clim = defaults[param]['clim']
    if not cmap:
        cmap = defaults[param]['cmap']

    # for electric field, multiply values by 1000 to get units of mV/m instead of V/m
    if param=='E':
        vmag = vmag*1000.
        dvmag = dvmag*1000.

    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(5,1)
    gs.update(hspace=0.4,left=0.15,right=0.9,bottom=0.05,top=0.95)

    for i, A in enumerate([vmag, dvmag, vdir, dvdir]):
        ax = plt.subplot(gs[i])
        f = ax.pcolormesh(A.T, vmin=clim[i][0], vmax=clim[i][1], cmap=plt.get_cmap(cmap[i]))
        ax.set_xticks(xticks)
        ax.set_xticklabels(time)
        ax.set_yticks(yticks)
        ax.set_yticklabels(binloc)
        ax.set_title(titles[i])
        pos = ax.get_position()
        cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
        cbar = fig.colorbar(f, cax=cax)

    ax = plt.subplot(gs[4])
    f = ax.quiver(vmag.T*np.sin(vdir.T*np.pi/180.), vmag.T*np.cos(vdir.T*np.pi/180.), np.sin(vdir.T*np.pi/180.), cmap=plt.get_cmap('coolwarm'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(np.arange(len(mlat)))
    ax.set_yticklabels(binloc)
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax, ticks=[-0.9,0,0.9])
    cbar.ax.set_yticklabels(['West','','East'])
    cbar.ax.quiverkey(f, 1.05, 1.1, clim[0][1], str(clim[0][1]))

    plt.savefig('{}magnitude.png'.format(param))
    plt.show()


def get_time_ticks(utime):
    # get x-axis (time) tick locations and labels
    # calculate datetime objects from unix time array
    time = [dt.datetime.utcfromtimestamp(np.mean(t)) for t in utime]
    # calculate experiment length
    exp_len = (time[-1] - time[0]).total_seconds()
    # set spaceing between each time tick based on lenght of experiment
    if exp_len < 1.*60.*60.:
        dtick = 10.
    elif exp_len < 3.*60.*60.:
        dtick = 15.
    elif exp_len < 6.*60.*60.:
        dtick = 30.
    else:
        dtick = 60.

    # get first time tick
    h0 = time[0].replace(minute=0, second=0)
    while h0<time[0]:
        h0 = h0+dt.timedelta(minutes=dtick)

    # get list of all time ticks
    time_ticks = []
    while h0 < time[-1]:
        time_ticks.append(h0)
        h0 = h0+dt.timedelta(minutes=dtick)

    # calculate position of each time tick
    xticks = [(tt-time[0]).total_seconds()/exp_len*len(time) for tt in time_ticks]
    # form list of strings for actual labels
    time = [t.strftime('%H:%M') for t in time_ticks]

    return time, xticks
