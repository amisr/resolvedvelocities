#!/usr/bin/env python

"""

"""

import os
import sys
import numpy as np
import tables
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib import cm
import matplotlib.dates

from datetime import datetime


class ResolvedAltPlotter(object):
    
    def __init__(self,file_to_plot):
        with tables.open_file(file_to_plot,'r') as h5:
            self.vest = h5.root.VectorVels.Vest.read()
            self.errvest = h5.root.VectorVels.errVest.read()
            self.alts = h5.root.VectorVels.Altitudes.read()
            self.times = h5.root.Time.UnixTime.read()
            self.inputfile = file_to_plot


    def make_plot(self,ylim=None,clim=1000.0,max_time=24.0,sc=1):

        cmap='jet'
        cmap_to_use = cm.get_cmap(cmap)
        cmap_to_use.set_bad('w',0)

        vmag = np.squeeze(np.sqrt(np.sum(self.vest**2,axis=2)))
        errvmag = np.squeeze(np.sqrt(np.sum(self.errvest**2,axis=2)))

        total_time = (self.times[-1,-1] - self.times[0,0]) / 3600.0
        num_time_groups = np.ceil(total_time / max_time)

        xlim = [self.times[0,0],self.times[-1,-1]]
        num_hours       = (xlim[-1] - xlim[0]) / 3600.0
        num_half_days   = num_hours / 12.0

        if num_half_days > 0.5:
            interval    = int(np.ceil(num_hours / 12.0))
            locator     = matplotlib.dates.HourLocator(interval=interval)
            formatter   = matplotlib.dates.DateFormatter("%H")
        elif num_half_days < 0.5:
            interval    = int(np.ceil(num_hours * 60.0 / 5.0 / sc))
            locator     = matplotlib.dates.MinuteLocator(interval=interval)
            formatter   = matplotlib.dates.DateFormatter("%H:%M")


        num_x, num_y    = vmag.shape
        temp_y          = self.alts
        temp_y          = np.repeat(temp_y[np.newaxis,:],num_x,axis=0)
        temp_y_diff     = np.repeat(np.diff(temp_y[0,:])[np.newaxis,:],num_x,axis=0)
        y_diff          = np.zeros(temp_y.shape)
        y_diff[:,0:-1]  = temp_y_diff
        y_diff[:,-1]    = temp_y_diff[:,-1]

        # Construct the range array for plotting
        y_plot              = np.zeros((num_x+1,num_y+1))
        y_plot[0:-1,0:-1]   = temp_y - y_diff/2
        y_plot[0:-1,-1]     = temp_y[:,-1] + y_diff[:,-1]/2
        y_plot[-1,:]        = y_plot[-2,:]

        if ylim is None:
            ylim = [y_plot[0,0],y_plot[0,-1]]

        # Construct the time array for plotting
        temp = np.zeros((num_x+1,num_y+1))
        temp[:num_x,:] = np.repeat(self.times[:,0][:,np.newaxis],num_y+1,axis=1)
        temp[num_x,:] = self.times[-1,-1]

        x_plot = np.array([[datetime.utcfromtimestamp(y) for y in x] for x in temp])
        xlim = [x_plot[0,0],x_plot[-1,0]]

        start_ind=0;
        for time_ind in range(int(num_time_groups)):
            end_ind = np.where(self.times[:,-1] <= (self.times[start_ind,0] + max_time * 3600.0))[0]
            end_ind = end_ind[-1]
            tlim = [start_ind,end_ind]
    
            if num_time_groups>1:
                txtra='_day' + str(time_ind)
            else:
                txtra=''
    
            # Figure out the time text for the title of the plots
            title =  "%s UT - "    % (xlim[0].strftime('%Y-%m-%d %H:%M:%S'))
            title += "%s UT"    % (xlim[1].strftime('%Y-%m-%d %H:%M:%S'))


            fig = pyplot.figure(figsize=(10,6))
            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            pc1 = ax1.pcolormesh(x_plot,y_plot,vmag,vmin=0,vmax=clim,cmap=cmap_to_use)
            pc2 = ax2.pcolormesh(x_plot,y_plot,errvmag,vmin=0,vmax=clim/10.,cmap=cmap_to_use)

            ax1.set_title(title)
            ax1.xaxis.set_major_locator(locator)
            ax1.xaxis.set_major_formatter(formatter)
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            # ax1.set_xticklabels([])
            ax2.xaxis.set_major_locator(locator)
            ax2.xaxis.set_major_formatter(formatter)
            ax2.set_xlim(xlim)
            ax2.set_ylim(ylim)
            # t = ax2.get_xticklabels()
            # ax2.set_xticklabels(t,rotation=90)

            ax1.set_ylabel('Altitude (km)')
            ax2.set_ylabel('Altitude (km)')
            ax2.set_xlabel('Time (UT)')

            cl1 = pyplot.colorbar(pc1,ax=ax1)
            cl1.set_label('Velocity Magnitude (m/s)')
            cl2 = pyplot.colorbar(pc2,ax=ax2)
            cl2.set_label('Velocity Error (m/s)')


            # Finally, save the figure if a name was provided.
            # Form the save name for the figure
            infilename = os.path.basename(self.inputfile)
            directory = os.path.dirname(self.inputfile)
            filename_wo_ext = os.path.splitext(infilename)[0]
            output_name = os.path.join(directory,filename_wo_ext+txtra+'.png')
            fig.savefig(output_name)


            start_ind = end_ind+1
            if start_ind >= self.times.shape[0]:
                break


if __name__ == '__main__':
    
    rap = ResolvedAltPlotter(sys.argv[1])
    rap.make_plot()