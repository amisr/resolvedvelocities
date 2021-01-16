#!/usr/bin/env python

"""

"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from matplotlib.colors import LogNorm

CMAP='RdBu'
rCMAP='RdBu_r'


# pad time gaps with NaNs for better pcolormesh plotting
def timegaps(time,data,rngOpt=[]):
    
    if len(rngOpt)>0:
        doRng=1
        rng2=[]
        
    time2=[]
    if np.ndim(data)==3:
        concnan=np.zeros((1,data.shape[1],data.shape[2]),dtype=data.dtype)*np.nan
    elif np.ndim(data)==2:
        concnan=np.zeros((1,data.shape[1]),dtype=data.dtype)*np.nan
    data2=data.copy()
    dt=np.median(np.diff(np.mean(time,axis=1)))
    for aa in range(time.shape[0]-1):
        time2.append(time[aa,0])
        if ( (time[aa+1,1]-time[aa,0]) > (dt*2.0) ):
            time2.append(time[aa,1])
            #print datetime.datetime.utcfromtimestamp(time[aa,1])
            if np.ndim(data)==3:
                data2=np.concatenate((data2[0:len(time2)-1,:,:],concnan,data2[len(time2)-1:,:,:]),axis=0)
            elif np.ndim(data)==2:
                data2=np.concatenate((data2[0:len(time2)-1,:],concnan,data2[len(time2)-1:,:]),axis=0)
    
    time2.append(time[-1,0])
    time2.append(time[-1,1])
    
    return np.array(time2), data2



class vvelsPlot:
    
    def __init__(self):
        """ initialization function """
        return
        
    def close(self):
        pyplot.close(self.figg)
        return
    
    def makePlot(self,time,MLTtime,lat,vx,vy,dvx,dvy,title='Vector Vels',units='m/s',parm='V',\
        p=[200.0,.25,4000,500.0],sc=15.0,cax=[-1000,1000],label='Mag. Lat. (degrees)',\
        ncols=2,vz=[],dvz=[],vzsc=1.0,nrows=4,geo=0,doQuiv=1,textsize=8,labsize=10):

        ncols=3
        nrows=3

        vx=vx.copy()
        vy=vy.copy()
        dvx=dvx.copy()
        dvy=dvy.copy()

        vz=vz.copy()*vzsc
        dvz=dvz.copy()*vzsc
        
        if lat.ndim==2:
            lat2=np.nanmean(lat,axis=1)
            lat=np.concatenate((lat[:,0],lat[[-1],1]))
        else:
            lat2=lat
        
        time2=np.mean(time,axis=1)

        textsize = 8        # size for axes text
        labsize = 10
        figBG   = 'w'        # the figure background color
        axesBG  = '#f6f6f6'  # the axies background color
        figsz = (10,9)


        figg=pyplot.figure(figsize=figsz, facecolor=figBG)
        ax = list()

        pos11 = [0.1  ,0.685,0.255,0.225]
        pos21 = [0.37 ,0.685,0.255,0.225]
        pos31 = [0.64 ,0.685,0.255,0.225]

        pos12 = [0.1  ,0.395,0.255,0.225]
        pos22 = [0.37 ,0.395,0.255,0.225]
        pos32 = [0.64 ,0.395,0.255,0.225]

        pos13 = [0.1  ,0.07 ,0.795,0.225]

        pos_col1 = [0.91  ,0.685,0.02,0.225]
        pos_col2 = [0.91  ,0.395 ,0.02,0.225]
        pos_col3 = [0.91  ,0.07 ,0.02,0.225]

        ax.append(figg.add_axes(pos11, facecolor=axesBG))
        ax.append(figg.add_axes(pos21, facecolor=axesBG))
        ax.append(figg.add_axes(pos31, facecolor=axesBG))

        ax.append(figg.add_axes(pos12, facecolor=axesBG))
        ax.append(figg.add_axes(pos22, facecolor=axesBG))
        ax.append(figg.add_axes(pos32, facecolor=axesBG))

        ax.append(figg.add_axes(pos13, facecolor=axesBG))

        ax.append(figg.add_axes(pos_col1, facecolor=axesBG))
        ax.append(figg.add_axes(pos_col2, facecolor=axesBG))
        ax.append(figg.add_axes(pos_col3, facecolor=axesBG))


        ii=0
        x,dat=timegaps(time,vx)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ylim=[np.nanmin(lat),np.nanmax(lat)]
        pc=ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=cax[0],vmax=cax[1],cmap=pyplot.get_cmap(CMAP))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_ylabel(label, fontsize=labsize)
        if geo==1:
            ax[ii].set_title('%s east (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        else:
            ax[ii].set_title('%s perp east (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        ax[ii].text(xlim[0],(ylim[1]-ylim[0])*0.15+ylim[1],title,fontsize=labsize, horizontalalignment='left')
        if ncols==1:
            cl=pyplot.colorbar(pc)
            cl.set_label(units,fontsize=labsize)
        

        ii=ii+1
        x,dat=timegaps(time,vy)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=cax[0],vmax=cax[1],cmap=pyplot.get_cmap(CMAP))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_yticklabels([], fontsize=textsize)
        if geo==1:
            ax[ii].set_title('%s north (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        else:
            ax[ii].set_title('%s perp north (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')

        ii=ii+1
        x,dat=timegaps(time,vz)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=cax[0],vmax=cax[1],cmap=pyplot.get_cmap(CMAP))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_yticklabels([], fontsize=textsize)
        if geo==1:
            ax[ii].set_title('%s up (%s) x %d' % (parm,units,vzsc), fontsize=labsize, horizontalalignment='center')	
        else:
            ax[ii].set_title('%s anti par (%s) x %d' % (parm,units,vzsc), fontsize=labsize, horizontalalignment='center')	
        
        ii=ii+1
        x,dat=timegaps(time,dvx)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        pc2=ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=0,vmax=cax[1]/5)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlabel('Time (UT)', fontsize=labsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_ylabel(label, fontsize=labsize)
        if geo==1:
            ax[ii].set_title('err %s east (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        else:
            ax[ii].set_title('err %s perp east (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')

        ii=ii+1
        x,dat=timegaps(time,dvy)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=0,vmax=cax[1]/5)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlabel('Time (UT)', fontsize=labsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_yticklabels([], fontsize=textsize)
        if geo==1:
            ax[ii].set_title('err %s north (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        else:
            ax[ii].set_title('err %s perp north (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
                
        ii=ii+1
        x,dat=timegaps(time,dvz)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=0,vmax=cax[1]/5)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlabel('Time (UT)', fontsize=labsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_yticklabels([], fontsize=textsize)
        if geo==1:
            ax[ii].set_title('err %s up (%s) x %d' % (parm,units,vzsc), fontsize=labsize, horizontalalignment='center')					
        else:
            ax[ii].set_title('err %s anti par (%s) x %d' % (parm,units,vzsc), fontsize=labsize, horizontalalignment='center')
        
        # quiver plot
        ii=ii+1
        I=np.where(np.isnan(vx))
        vx[I]=0
        vy[I]=0
        I=np.where(np.isnan(vy))
        vx[I]=0
        vy[I]=0
        
        I=np.where(np.absolute(dvx)/(np.absolute(vx)+p[0])>p[1])
        vx[I]=0
        vy[I]=0
        I=np.where(np.absolute(dvy)/(np.absolute(vy)+p[0])>p[1])
        vx[I]=0
        vy[I]=0
        

        I=np.where(np.sqrt(vx*vx+vy*vy)>p[2])
        vx[I]=0
        vy[I]=0
    
        I=np.where((np.absolute(dvx)>p[3]) | (np.absolute(dvy)>p[3]))
        vx[I]=0
        vy[I]=0
        
        C=np.ones(np.transpose(vx).shape)-np.sign(np.transpose(vx))
        x=matplotlib.dates.epoch2num(time2)
        [X,Y]=np.meshgrid(x,lat2)
        Q=ax[ii].quiver(X,Y,np.transpose(vx),np.transpose(vy),C,clim=[0,2],scale=900.0*sc,width=0.3*0.005,cmap=pyplot.get_cmap(rCMAP))
        ax[ii].quiverkey(Q, 1.06, 1.15, cax[1], str(cax[1]) + ' ' + units,fontproperties={'size' : labsize},labelpos='S')
        ax[ii].xaxis.tick_bottom()
        ax[ii].set_xlabel('Time (UT)', fontsize=labsize)
        labels = ax[ii].get_xticklabels()
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_ylabel(label, fontsize=labsize)
        ax[ii].set_xlim((x[0],x[-1]))
        ax[ii].set_ylim((np.nanmin(lat2),np.nanmax(lat2)))	
        
        ax22 = figg.add_axes(ax[ii].get_position(), sharey=ax[ii], frameon=False)
        ax22.plot(MLTtime,np.nan*np.array(MLTtime),'k')
        ax22.set_xlim([MLTtime[0],MLTtime[-1]])
        ax22.xaxis.tick_top()
        ax22.set_xlabel('Magnetic Local Time', fontsize=labsize)
        ax22.xaxis.set_label_position('top')
        ax22.tick_params(axis='both',labelsize=textsize)
        ax22.set_ylim((np.nanmin(lat2),np.nanmax(lat2)))

        # Colorbar top
        ii=ii+1
        cl=pyplot.colorbar(pc,ax[ii])
        cl.set_label(units,fontsize=labsize)
        labels = ax[ii].get_yticklabels()
        ax[ii].set_yticklabels(labels, fontsize=textsize)

        # Colorbar middle
        ii=ii+1
        cl=pyplot.colorbar(pc2,ax[ii])
        cl.set_label(units,fontsize=labsize)
        labels = ax[ii].get_yticklabels()
        ax[ii].set_yticklabels(labels, fontsize=textsize)

        # Colorbar bottom
        ii=ii+1
        cl=pyplot.colorbar(Q,ax[ii])
        ax[ii].invert_yaxis()
        cl.set_ticks(np.linspace(0,2, 5))
        cl.set_ticklabels(['','East', '', 'West',''])

        # Set time formatting
        locator = matplotlib.dates.HourLocator(interval=1)
        formatter = matplotlib.dates.DateFormatter("%H:%M")
        
        dx=(time2[-1]-time2[0])/3600.0
        dx2=dx/7.0
        

        ss=7
        for rr in range(len(ax)):
        
            if dx2>0.5:
                interval=int(np.ceil(dx/7.0))
                locator = matplotlib.dates.HourLocator(interval=interval)
                formatter = matplotlib.dates.DateFormatter("%H:%M")
            elif dx2<0.5:
                interval=int(np.ceil(dx*60.0/7.0))
                locator = matplotlib.dates.MinuteLocator(interval=interval)
                formatter = matplotlib.dates.DateFormatter("%H:%M")
                
            ax[rr].xaxis.set_major_locator(locator)
            ax[rr].xaxis.set_major_formatter(formatter)

        self.figg=figg
        
        return
  
class vvelsMagPlot:
    
    def __init__(self):
        """ initialization function """
        return
        
    def close(self):
        pyplot.close(self.figg)
        return        
                    
    def makePlot(self,time,MLTtime,lat,vx,vy,dvx,dvy,title='Vector Vels',units='m/s',parm='V',cax1=[-1000,1000],cax2=[-180,180],label='Mag. Lat. (degrees)'):

        vx=vx.copy()
        vy=vy.copy()
        dvx=dvx.copy()
        dvy=dvy.copy()
        
        if lat.ndim==2:
            lat2=np.nanmean(lat,axis=1)
            lat=np.concatenate((lat[:,0],lat[[-1],1]))
        else:
            lat2=lat
        
        time2=np.mean(time,axis=1)

        textsize = 8        # size for axes text
        labsize = 10
       
        figBG   = 'w'        # the figure background color
        axesBG  = '#f6f6f6'  # the axies background color
        figsz = (7,9)

        dx, dy= 0.015, 0.05	
        nrows=4; ncols=1
        POS=[0.07,0.75,0.8-dx,1.0/(nrows)-dy*1.5]
        POS=[0.1,0.75,0.77-dx,1.0/(nrows)-dy*1.5]
            
        figg=pyplot.figure(figsize=figsz, facecolor=figBG)
        
        ax=[]
        for aa in range(nrows):
            for bb in range(ncols):
                rect=[POS[0]+(POS[2]+dx)*bb,POS[1]-(POS[3]+dy)*aa,POS[2],POS[3]]
                ax.append(figg.add_axes(rect, facecolor=axesBG))
        
        pc=[]
        
        ii=0
        x,dat=timegaps(time,vx)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        ylim=[np.nanmin(lat),np.nanmax(lat)]
        pc.append(ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=cax1[0],vmax=cax1[1]))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        # labels = ax[ii].get_xticklabels()
        # ax[ii].set_xticklabels(labels, fontsize=textsize)
        # labels = ax[ii].get_yticklabels()
        # ax[ii].set_yticklabels(labels, fontsize=textsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)	
        ax[ii].set_ylabel(label, fontsize=labsize)
        ax[ii].set_title('%s mag (%s)' % (parm,units), fontsize=labsize, horizontalalignment='center')
        ax[ii].text(xlim[0],(ylim[1]-ylim[0])*0.15+ylim[1],title,fontsize=labsize, horizontalalignment='left')
                    
        ii=ii+1
        
        x,dat=timegaps(time,dvx)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        pc.append(ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=0.0,vmax=cax1[1]/5))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        # labels = ax[ii].get_xticklabels()
        # ax[ii].set_xticklabels(labels, fontsize=textsize)
        # labels = ax[ii].get_yticklabels()
        # ax[ii].set_yticklabels(labels, fontsize=textsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)
        ax[ii].set_ylabel(label, fontsize=labsize)
        ax[ii].set_title('err %s mag (%s)' %(parm,units), fontsize=labsize, horizontalalignment='center')

        ii=ii+1
        
        x,dat=timegaps(time,vy)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        pc.append(ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=cax2[0],vmax=cax2[1],cmap=pyplot.get_cmap('hsv')))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        # labels = ax[ii].get_xticklabels()
        # ax[ii].set_xticklabels(labels, fontsize=textsize)
        # labels = ax[ii].get_yticklabels()
        # ax[ii].set_yticklabels(labels, fontsize=textsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)	
        ax[ii].set_ylabel(label, fontsize=labsize)
        ax[ii].set_title('%s dir (%s)' % (parm,'degrees'), fontsize=labsize, horizontalalignment='center')

        ii=ii+1
        
        x,dat=timegaps(time,dvy)
        dat=np.ma.masked_where(np.isnan(dat),dat)
        x=matplotlib.dates.epoch2num(x)
        xlim=[x[0],x[-1]]	
        pc.append(ax[ii].pcolor(x,lat,np.transpose(dat),edgecolors='none',vmin=0,vmax=cax2[1]/5))
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlabel('Time (UT)', fontsize=labsize)
        # labels = ax[ii].get_xticklabels()
        # ax[ii].set_xticklabels(labels, fontsize=textsize)
        # labels = ax[ii].get_yticklabels()
        # ax[ii].set_yticklabels(labels, fontsize=textsize)
        ax[ii].tick_params(axis='both',labelsize=textsize)	
        ax[ii].set_ylabel(label, fontsize=labsize)
        ax[ii].set_title('err %s dir (%s)' % (parm,'degrees'), fontsize=labsize, horizontalalignment='center')
                
        ii=ii+1
        
        if ncols==1:
            bb=1
            for aa in range(len(pc)):
                rect=[POS[0]+(POS[2]+dx)*bb,POS[1]-(POS[3]+dy)*aa,POS[2]/20,POS[3]]
                ax.append(figg.add_axes(rect, facecolor=axesBG))	
                cl=pyplot.colorbar(pc[aa],ax[ii])
                if aa<2:
                    cl.set_label(units,fontsize=labsize)
                else:
                    cl.set_label('degrees',fontsize=labsize)
                labels = ax[ii].get_yticklabels()
                ax[ii].set_yticklabels(labels, fontsize=textsize)
                ii+=1
                    
        locator = matplotlib.dates.HourLocator(interval=1)
    #	locator = matplotlib.dates.MinuteLocator(interval=5)
        formatter = matplotlib.dates.DateFormatter("%H:%M")
        
        dx=(time2[-1]-time2[0])/3600.0
        dx2=dx/7.0
        
        for rr in range(len(ax)):
        
            if dx2>0.5:
                interval=int(np.ceil(dx/7.0))
                locator = matplotlib.dates.HourLocator(interval=interval)
                formatter = matplotlib.dates.DateFormatter("%H:%M")
            elif dx2<0.5:
                interval=int(np.ceil(dx*60.0/7.0))
                locator = matplotlib.dates.MinuteLocator(interval=interval)
                formatter = matplotlib.dates.DateFormatter("%H:%M")
                
            ax[rr].xaxis.set_major_locator(locator)
            ax[rr].xaxis.set_major_formatter(formatter)
        
        # pyplot.show()
        
        self.figg=figg
        
        return figg
