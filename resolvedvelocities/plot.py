# plot.py

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import datetime as dt

def plot_components(utime, mlat, mlon, vector, covariance, param='V'):

    # get x-axis (time) tick locations and labels
    time, xticks = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    binloc = ['({} N,\n {} E)'.format(lat, lon) for lat, lon in zip(mlat, mlon)]

    if param=='V':
        # get plot titles
        vtitles = ['Ve1 (m/s)', 'Ve2 (m/s)', 'Ve3 (m/s) x 10']
        etitles = ['dVe1 (m/s)', 'dVe2 (m/s)', 'dVe3 (m/s) x 10']
        # get plot limits
        vlim = [-1500., 1500.]
        elim = [0., 350.]
    elif param=='E':
        # get plot titles
        vtitles = ['Ed1 (mV/m)', 'Ed2 (mV/m)', 'Ed3 (mV/m) x 10']
        etitles = ['dEd1 (mV/m)', 'dEd2 (mV/m)', 'dEd3 (mV/m) x 10']
        # get plot limits
        vlim = [-80., 80.]
        elim = [0., 15.]
        # multiply values by 1000. to get units of mV/m instead of V/m
        vector = vector*1000.
        covariance = covariance*1000.*1000.

    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(2,3)
    gs.update(wspace=0.3,hspace=0.2,left=0.05,right=0.9,bottom=0.05,top=0.95)

    # TODO: These kinds of plots don't really make much sense if we're allowing arbitrary bin locations
    # make figure save, not show to screen
    for i in range(3):

        # for along B component, multiply by 10
        if i==2:
            fac=10.
        else:
            fac=1.

        ax1 = plt.subplot(gs[0,i])
        f1 = ax1.pcolormesh(vector[:,:,i].T*fac, vmin=vlim[0], vmax=vlim[1], cmap=plt.get_cmap('coolwarm'))
        ax1.set_xticks(xticks)
        ax1.set_xticklabels(time)
        ax1.set_yticks(yticks)
        ax1.set_yticklabels(binloc)
        ax1.set_title(vtitles[i])
        ax1.tick_params(labelsize=8)

        ax2 = plt.subplot(gs[1,i])
        f2 = ax2.pcolormesh(np.sqrt(covariance[:,:,i,i]).T*fac, vmin=elim[0], vmax=elim[1])
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(time)
        ax2.set_yticks(yticks)
        ax2.set_yticklabels(binloc)
        ax2.set_xlabel('UT')
        ax2.set_title(etitles[i])
        ax2.tick_params(labelsize=8)

    pos = ax1.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f1, cax=cax)

    pos = ax2.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f2, cax=cax)

    plt.savefig('{}components.png'.format(param))
    plt.show()


def plot_magnitude(utime, mlat, mlon, vmag, dvmag, vdir, dvdir, param='V'):

    # get x-axis (time) tick locations and labels
    time, xticks = get_time_ticks(utime)

    # get y-axis (bin) tick locations and labels
    yticks = np.arange(len(mlat))+0.5
    binloc = ['({} N, {} E)'.format(lat, lon) for lat, lon in zip(mlat, mlon)]

    # set titles
    if param=='V':
        titles = ['V magnitude (m/s)', 'V magnitude error (m/s)', 'V direction (deg)', 'V direction error (deg)', '']
        mlim = [0.,1500.]
        emlim = [0., 350.]
    elif param=='E':
        titles = ['E magnitude (mV/m)', 'V magnitude error (mV/m)', 'E direction (deg)', 'E direction error (deg)', '']
        vmag = vmag*1000.
        dvmag = dvmag*1000.
        mlim = [0.,80.]
        emlim = [0., 15.]


    dlim = [-180., 180.]
    edlim = [0., 35.]

    # TODO: consolodate this into a loop of some kind
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(5,1)
    gs.update(hspace=0.4,left=0.15,right=0.9,bottom=0.05,top=0.95)

    ax = plt.subplot(gs[0])
    f = ax.pcolormesh(vmag.T, vmin=mlim[0], vmax=mlim[1])
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(yticks)
    ax.set_yticklabels(binloc)
    ax.set_title(titles[0])
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax)

    ax = plt.subplot(gs[1])
    f = ax.pcolormesh(dvmag.T, vmin=emlim[0], vmax=emlim[1])
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(yticks)
    ax.set_yticklabels(binloc)
    ax.set_title(titles[1])
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax)

    ax = plt.subplot(gs[2])
    f = ax.pcolormesh(vdir.T, vmin=dlim[0], vmax=dlim[1], cmap=plt.get_cmap('hsv'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(yticks)
    ax.set_yticklabels(binloc)
    ax.set_title(titles[2])
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax)

    ax = plt.subplot(gs[3])
    f = ax.pcolormesh(dvdir.T, vmin=edlim[0], vmax=edlim[1])
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(yticks)
    ax.set_title(titles[3])
    ax.set_yticklabels(binloc)
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax)

    ax = plt.subplot(gs[4])
    f = ax.quiver(vmag.T*np.sin(vdir.T*np.pi/180.), vmag.T*np.cos(vdir.T*np.pi/180.), np.sin(vdir.T*np.pi/180.), cmap=plt.get_cmap('coolwarm'))
    ax.set_xticks(xticks)
    ax.set_xticklabels(time)
    ax.set_yticks(np.arange(len(mlat)))
    ax.set_title(titles[4])
    ax.set_yticklabels(binloc)
    pos = ax.get_position()
    cax = fig.add_axes([0.91, pos.y0, 0.015, pos.y1-pos.y0])
    cbar = fig.colorbar(f, cax=cax, ticks=[-0.9,0,0.9])
    cbar.ax.set_yticklabels(['West','','East'])
    cbar.ax.quiverkey(f, 1.05, 1.1, mlim[1], str(mlim[1]))

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
