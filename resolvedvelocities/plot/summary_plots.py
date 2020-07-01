# summary_plots.py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import datetime as dt
import os

def plot_components(utime, mlat, mlon, vector, covariance, param='V', titles=None, clim=None, cmap=None, savedir=None):

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
    if not savedir:
        savedir = os.getcwd()

    # for electric field, multiply values by 1000 to get units of mV/m instead of V/m
    if param=='E':
        vector = vector*1000.
        covariance = covariance*1000.*1000.

    # pad time gaps
    utime, [vector, covariance] = timegaps(utime, [vector, covariance])

    # get x-axis (time) tick locations and labels
    time, xticks, xlims = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    binmlat = ['{:.2f} N'.format(ml) for ml in mlat]

    xedge = np.append(utime[:,0],utime[-1,1])
    yedge = np.arange(len(mlat)+1)

    # loop over every day of experiment
    for xlim in xlims:

        fig = plt.figure(figsize=(15,10))
        gs = gridspec.GridSpec(2,3)
        gs.update(wspace=0.3,hspace=0.2,left=0.05,right=0.9,bottom=0.05,top=0.9)

        for i, A in enumerate([vector,np.sqrt(np.diagonal(covariance,axis1=-1,axis2=-2))]):
            for j in range(3):
                # for along B component, multiply by 10
                if j==2:
                    fac=10.
                else:
                    fac=1.

                ax = plt.subplot(gs[i,j])
                f = ax.pcolormesh(xedge, yedge, A[:,:,j].T*fac, vmin=clim[i][0], vmax=clim[i][1], cmap=plt.get_cmap(cmap[i]))
                ax.set_xticks(xticks)
                ax.set_xticklabels(time)
                ax.set_yticks(yticks)
                ax.set_yticklabels(binmlat)
                ax.set_xlim(xlim)
                ax.set_title(titles[i][j])
                ax.set_xlabel('Universal Time')
                ax.set_ylabel('Apex MLAT')
                ax.tick_params(labelsize=8)

            pos = ax.get_position()
            cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
            cbar = fig.colorbar(f, cax=cax)

        datestr = '{:%Y-%m-%d %H:%M:%S} - {:%Y-%m-%d %H:%M:%S}'.format(dt.datetime.utcfromtimestamp(xlim[0]),dt.datetime.utcfromtimestamp(xlim[1]))
        plt.gcf().text(0.02, 0.95, datestr, fontsize=12)

        filename = os.path.join(savedir,'{}compontents_{:%Y%m%d%H%M%S}_{:%Y%m%d%H%M%S}.png'.format(param,dt.datetime.utcfromtimestamp(xlim[0]),dt.datetime.utcfromtimestamp(xlim[1])))
        plt.savefig(filename)
        # plt.show()
        plt.close()


def plot_magnitude(utime, mlat, mlon, vmag, dvmag, vdir, dvdir, param='V', titles=None, clim=None, cmap=None, savedir=None):

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
    if not savedir:
        savedir = os.getcwd()

    # for electric field, multiply values by 1000 to get units of mV/m instead of V/m
    if param=='E':
        vmag = vmag*1000.
        dvmag = dvmag*1000.

    # pad time gaps
    utime, [vmag, dvmag, vdir, dvdir] = timegaps(utime, [vmag, dvmag, vdir, dvdir])

    # get x-axis (time) tick locations and labels
    time, xticks, xlims = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    # binloc = ['({:.1f} N, {:.1f} E)'.format(lat, lon) for lat, lon in zip(mlat, mlon)]
    binmlat = ['{:.2f} N'.format(ml) for ml in mlat]

    xedge = np.append(utime[:,0],utime[-1,1])
    yedge = np.arange(len(mlat)+1)

    # loop over every day of experiment
    for xlim in xlims:

        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(5,1)
        gs.update(hspace=0.4,left=0.15,right=0.9,bottom=0.05,top=0.94)

        for i, A in enumerate([vmag, dvmag, vdir, dvdir]):
            ax = plt.subplot(gs[i])
            f = ax.pcolormesh(xedge, yedge, A.T, vmin=clim[i][0], vmax=clim[i][1], cmap=plt.get_cmap(cmap[i]))
            ax.set_xticks(xticks)
            ax.set_xticklabels(time)
            ax.set_yticks(yticks[::2])
            ax.set_yticklabels(binmlat[::2])
            ax.set_xlim(xlim)
            # ax.set_xlabel('Universal Time')
            ax.set_ylabel('Apex MLAT')
            ax.set_title(titles[i])
            pos = ax.get_position()
            cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
            cbar = fig.colorbar(f, cax=cax)

        ax = plt.subplot(gs[4])
        f = ax.quiver(xedge[:-1], yedge[:-1], vmag.T*np.sin(vdir.T*np.pi/180.), vmag.T*np.cos(vdir.T*np.pi/180.), np.sin(vdir.T*np.pi/180.), cmap=plt.get_cmap('coolwarm'))
        ax.set_xticks(xticks)
        ax.set_xticklabels(time)
        ax.set_yticks(np.arange(len(mlat))[::2])
        ax.set_yticklabels(binmlat[::2])
        ax.set_xlim(xlim)
        ax.set_xlabel('Universal Time')
        ax.set_ylabel('Apex MLAT')
        pos = ax.get_position()
        cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
        cbar = fig.colorbar(f, cax=cax, ticks=[-0.9,0,0.9])
        cbar.ax.set_yticklabels(['West','','East'])
        cbar.ax.quiverkey(f, 1.05, 1.1, clim[0][1], str(clim[0][1]))

        datestr = '{:%Y-%m-%d %H:%M:%S} - {:%Y-%m-%d %H:%M:%S}'.format(dt.datetime.utcfromtimestamp(xlim[0]),dt.datetime.utcfromtimestamp(xlim[1]))
        plt.gcf().text(0.02, 0.97, datestr, fontsize=12)

        filename = os.path.join(savedir,'{}magnitude_{:%Y%m%d%H%M%S}_{:%Y%m%d%H%M%S}.png'.format(param,dt.datetime.utcfromtimestamp(xlim[0]),dt.datetime.utcfromtimestamp(xlim[1])))
        plt.savefig(filename)
        # plt.show()
        plt.close()


def get_time_ticks(utime):
    # get x-axis (time) tick locations and labels

    # calculate experiment length
    exp_start = dt.datetime.utcfromtimestamp(utime[0,0])
    exp_end = dt.datetime.utcfromtimestamp(utime[-1,1])
    exp_len = (exp_end-exp_start).total_seconds()
    # set spaceing between each time tick based on lenght of experiment
    if exp_len < 1.*60.*60.:
        dtick = 10.
    elif exp_len < 3.*60.*60.:
        dtick = 15.
    elif exp_len < 6.*60.*60.:
        dtick = 30.
    elif exp_len < 12.*60.*60.:
        dtick = 60.
    else:
        dtick = 3*60.

    # get first time tick
    h0 = exp_start.replace(minute=0, second=0)
    while h0<exp_start:
        h0 = h0+dt.timedelta(minutes=dtick)

    # get list of all time ticks
    time_ticks = []
    while h0 < exp_end:
        time_ticks.append(h0)
        h0 = h0+dt.timedelta(minutes=dtick)

    # calculate position of each time tick
    xticks = [(tt-dt.datetime.utcfromtimestamp(0)).total_seconds() for tt in time_ticks]
    # form list of strings for actual labels
    time = [t.strftime('%H:%M') for t in time_ticks]

    xlims = []
    if exp_len > 24.*60.*60:
    # if experiment is longer than 1 day, break up plot limits and create a series of plots
        st = exp_start
        et = exp_start.replace(hour=0, minute=0, second=0) + dt.timedelta(days=1)
        while et < exp_end:
            xlims.append([st,et])
            st = et
            et = et + dt.timedelta(days=1)
        xlims.append([st,exp_end])
    else:
    # if experiment shorter than one day, just create one plot
        xlims.append([exp_start, exp_end])

    # convert xlims to unix time
    xlims = [[(x[0]-dt.datetime.utcfromtimestamp(0)).total_seconds(), (x[1]-dt.datetime.utcfromtimestamp(0)).total_seconds()] for x in xlims]

    return time, xticks, xlims


def timegaps(time, data_arrays):
# pad time gaps with NaNs for better pcolormesh plotting

    # find the time between sucessive records and the cadance of the entire time series
    time_diff = np.diff(np.mean(time,axis=1))
    dt = np.median(time_diff)
    # find gaps in the time series
    gaps = np.argwhere(time_diff > 2*dt).flatten()+1

    # if no gaps, return original arrays
    if not gaps.size:
        return time, data_arrays

    # create array of times to fill each gap
    time_insert = np.array([time[g-1]+dt for g in gaps])
    # insert fill times into each gap
    time2 = np.insert(time,gaps,time_insert, axis=0)

    data_arrays2 = []
    for data in data_arrays:
        # create array of nans to fill each gap in data array
        data_insert = np.full(data[0].shape,np.nan)
        # insert fill nans into data array
        data2 = np.insert(data,gaps,data_insert, axis=0)
        data_arrays2.append(data2)

    return time2, data_arrays2
