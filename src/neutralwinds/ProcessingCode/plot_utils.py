#! /usr/bin/env python

"""
xxxxx

~M. Nicolls
last revised: xx/xx/2007

"""

import os
import scipy, tables
import matplotlib
matplotlib.use('Agg')
import pylab, matplotlib.dates

class MyClass:
    "A simple example class"
    i = 12345
    def f(self):
        return 'hello world'

def geoplot(az,el,rng,ht,plat,plong,dip,dec,Ibeam=[]):

    labsize=10
    textsize=8

    I=scipy.where(el>0)[0]
    az=az[I]
    el=el[I]
    plat=plat[I]
    plong=plong[I]
    dip=dip[I]
    dec=dec[I]
#   ht=ht[I,:]

    Nbeams=az.size

    elm=int(scipy.floor(el.min()/10.0)*10.0)

    axesBG  = '#f6f6f6'
    figg=pylab.figure(1, figsize=(5,8), facecolor='w')

    # polar plot
    r = (90.0-el)*pi/180
    theta = (-az+90)*pi/180
    area = 50
    colors = r
    ax = pylab.subplot(211, polar=True,axisbg=axesBG)
#    c = pylab.scatter(theta, r, c=colors, s=area)
    c = pylab.scatter(theta, r, s=area)
    c.set_alpha(0.75)

    ring_angles = [ x * (pi/180.) for x in range(0,100-elm,10)]
    ring_labels = [ str(x) for x in range(elm,100,10)]
    ring_labels.reverse()
    if ring_angles[0]==0.0:
        ring_angles[0]=ring_angles[0]+0.01
    pylab.rgrids (ring_angles,ring_labels, angle=-90+22.5,fontsize=textsize)
    pylab.thetagrids( range(0,360,45), ('E', 'NE', 'N','NW', 'W', 'SW', 'S', 'SE'), fontsize=textsize )

    if len(Ibeam)==0:
        Ibeam=range(Nbeams)

    for ii in range(az.size):
#       t=r'$%d: \ %2.1f^o,\ %2.1f^o$' % (ii,az[ii],el[ii])
        ri=scipy.where(scipy.array(Ibeam)==ii)[0]
        t=r'$%d$' % (ri+1)
        #print ri

        ax.text(theta[ii],r[ii],t,fontsize=labsize, horizontalalignment='right')


    #

    ax = pylab.subplot(212,axisbg=axesBG)
    for ii in range(Nbeams):
        ri=scipy.where(scipy.array(Ibeam)==ii)[0]
        t=r'$%d$' % (ri+1)
        tlat=plat[ii,:].copy()
        tlon=plong[ii,:].copy()
        I=scipy.where(scipy.absolute(scipy.diff(tlon))>300.0)[0]
        if len(I)>0:
            tlon[I[0]+1:]+=360.0
#       pylab.plot(plat[ii,:],plong[ii,:],label=str(ii),c=r[ii],linewidth=2)
#        pylab.plot(plat[ii,:],plong[ii,:],label=str(ii),c='%f' % r[ii],linewidth=2)
        pylab.plot(tlat,tlon,label=str(ii),linewidth=2)
        ax.text(tlat[-1],tlon[-1],t,fontsize=labsize, horizontalalignment='right')

    labels = pylab.getp(ax, 'xticklabels')
    pylab.setp(labels, fontsize=textsize)
    labels = pylab.getp(ax, 'yticklabels')
    pylab.setp(labels, fontsize=textsize)

    ax.set_ylabel('Mag. Long. (degrees)', fontsize=labsize)
    ax.set_xlabel('Mag. Lat. (degrees)', fontsize=labsize)


    ##pylab.show()

    return figg

def multi_axes(nrows,ncols,dx=0.02,dy=0.05,figsz=(14,10)):

    pylab.ioff()

    figBG   = 'w'        # the figure background color
    axesBG  = '#f6f6f6'  # the axies background color

    if ncols==1 and nrows==1:
        POS=[0.1,0.1,1.0/(ncols+3.0)-dx,1.0/(nrows)-dy*5]
    elif ncols==1:
        POS=[0.1,0.65,1.0/(ncols+3.0)-dx,1.0/(nrows)-dy*2]
    elif nrows==1:
        POS=[0.1,0.1,1.0/(ncols+1)-dx,1.0/(nrows)-dy*5]
    elif ncols==-1 and nrows==-1:
        ncols=1; nrows=2; dx=0.1
        POS=[0.1,0.1,1.0-dx,1.0/(nrows)-dy*5]
    else:
        fy=1.0/(nrows)-dy*1.5
        POS=[0.1,1.0-fy-dy*1.5,1.0/(ncols+1)-dx,fy]

    figg=pylab.figure(figsize=figsz, facecolor=figBG)

    ax=[]
    for aa in range(nrows):
        for bb in range(ncols):
            rect=[POS[0]+(POS[2]+dx)*bb,POS[1]-(POS[3]+dy)*aa,POS[2],POS[3]]
            ax.append(pylab.axes(rect))

    return figg,ax

def plot_array(ax,nrows,ncols,*args):

    # arr should be

    pylab.ioff()

    textsize = 8       # size for axes text
    labsize = 10

    xlims=args[-5]
    xlabels=args[-4]
    ylabels=args[-3]
    title=args[-2]
    bmcodes=args[-1]

    ii=0
    for aa in range(nrows):
        I=scipy.where(scipy.isfinite(args[aa*3+1]))

        ymin=args[aa*3+1][I].min()-15
        ymax=args[aa*3+1][I].max()+15

        for bb in range(ncols):

            ax[ii].clear()
            ax[ii].hold(True)
            if aa==0:
                ax[ii].set_xscale("log")
            for cc in range(args[aa*3].shape[2]):
                tval=args[aa*3][bb,:,cc]
                tax=args[aa*3+1][bb,:,cc]
                terr=scipy.repeat(args[aa*3+2][bb,:,cc][scipy.newaxis,:],2,axis=0)
                if aa==0:
                    I=scipy.where(tval<=0.0)[0]
                    tval[I]=scipy.nan
                    terr[:,I]=scipy.nan
                    I=scipy.where(terr[0,:]>=tval)[0]
                    terr[0,I]=tval[I]-10.0**-10
                I=scipy.where((tval<xlims[aa][0]) | (tval>xlims[aa][1]))[0]
                tval[I]=scipy.nan; terr[:,I]=scipy.nan
                I=scipy.where((tval-terr[0,:])<xlims[aa][0])[0]; terr[0,I]=tval[I]-xlims[aa][0]-1.0e-10
                I=scipy.where((tval+terr[1,:])>xlims[aa][1])[0]; terr[1,I]=xlims[aa][1]-tval[I]-1.0e-10
                ax[ii].errorbar(tval,tax,None,terr,'-k.')
            if aa==0 and bb==0:
                ax[ii].text(1,(ymax-ymin)*0.15+ymax,title,fontsize=labsize, horizontalalignment='left')
            if bb>0:
                ax[ii].yaxis.set_ticklabels([])
            else:
                ax[ii].set_ylabel(ylabels[aa], fontsize=labsize)
                labels = pylab.getp(ax[ii], 'yticklabels')
                pylab.setp(labels, fontsize=textsize)

            ax[ii].set_xlabel(xlabels[aa], fontsize=labsize)
            ax[ii].set_xlim(xlims[aa])
            ax[ii].set_ylim((ymin,ymax))

            labels = pylab.getp(ax[ii], 'xticklabels')
            pylab.setp(labels, fontsize=textsize)
            if aa==0:
                tt=r"$(%.1f^o \rm{az} , \ %.1f^o \rm{el})$" % (bmcodes[ii,1],bmcodes[ii,2])
                ax[ii].set_title(tt, fontsize=labsize, horizontalalignment='center')
            ii=ii+1

    return

def pcolor_plot_dual(x,y,y2,data,cax,xlim,ylim,xl,yl,yl2,title,text,bmcodes,log=0,ax=[],figg=-543,xxx=0.0,yyy=0.0,sc=1.0,im=0):

    pylab.ioff()

    textsize = 8        # size for axes text
    labsize = 12

    Nbeams=data.shape[1]
    Nplots=Nbeams+1
    if Nbeams==1:
        nrows=1
        ncols=2
    else:
        nrows=int(scipy.ceil(scipy.sqrt(Nplots)))
        ncols=nrows
        if nrows*(nrows-1)>=Nplots:
            ncols=nrows-1

    dx=(x[-1]-x[0])/3600.0
    dx2=dx/12.0
    if dx2>0.5:
        interval=scipy.ceil(dx/12.0)
        locator = matplotlib.dates.HourLocator(interval=interval)
        formatter = matplotlib.dates.DateFormatter("%H")
    elif dx2<0.5:
        interval=scipy.ceil(dx*60.0/3.0/sc)
        locator = matplotlib.dates.MinuteLocator(interval=interval)
        formatter = matplotlib.dates.DateFormatter("%H:%M")

    x=matplotlib.dates.epoch2num(x)
    xlim=[x[0],x[-1]]

    if len(ax)==0 or pylab.gcf()!=figg:
        (figg,ax)=multi_axes(nrows,ncols)

    for ii in range(Nbeams):
        ax[ii].clear()
        pc=ax[ii].pcolor(x,y[ii,:],scipy.transpose(data[:,ii,:]),vmin=cax[0],vmax=cax[1])
        ax[ii].xaxis.set_major_locator(locator)
        ax[ii].xaxis.set_major_formatter(formatter)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        if scipy.mod(ii,ncols)==0:
            ax[ii].set_ylabel(yl, fontsize=labsize)
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        labels = pylab.getp(ax[ii], 'yticklabels')
        pylab.setp(labels, fontsize=textsize)
#       if ii>=(nrows*(ncols-1)):
        if scipy.mod(ii,nrows)==0:
            ax[ii].set_xlabel(xl, fontsize=labsize)
        tt=r"$(%.1f^o \ \rm{az} , \ %.1f^o \ \rm{el})$" % (bmcodes[ii,1],bmcodes[ii,2])
        ax[ii].set_title(tt, fontsize=labsize, horizontalalignment='center')
        if ii==0:
            ax[ii].text(xlim[0]+(xlim[1]-xlim[0])/2+xxx,(ylim[1]-ylim[0])*0.05*nrows+ylim[1]+yyy,title,fontsize=labsize, horizontalalignment='left')

    if len(ax)==Nplots:
        ii=ii+1
        try:
            tmp=ax[ii].get_position()
            gp=scipy.array([tmp.xmin,tmp.ymin,tmp.width,tmp.height])
        except:
            gp=ax[ii].get_position()
    #   gp[1]=gp[1]*2.
        gp[3]=gp[3]/6.
        ax[ii].set_position(gp)
        if log:
            cl=pylab.colorbar(pc,ax[ii],orientation='horizontal',format=pylab.FormatStrFormatter('$10^{%.1f}$'))
        else:
            cl=pylab.colorbar(pc,ax[ii],orientation='horizontal')
        ax[ii].yaxis.set_ticklabels([])
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
    #   labels=scipy.linspace(cax[0],cax[1],len(labels)).tolist()
        cl.set_label(text,fontsize=labsize)
    else:
        if log:
            cl=figg.colorbar(pc,orientation='vertical',format=pylab.FormatStrFormatter('$10^{%.1f}$'),pad=0.15)
        else:
            cl=figg.colorbar(pc,orientation='vertical',pad=0.15)
        cl.set_label(text,fontsize=labsize*1.25)

    for jj in range(len(ax)):
        ax2=pylab.twinx(ax[ii])
        labels = pylab.getp(ax2, 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
 #       ax2.plot([x[0],x[-1]],[y2[0],y2[-1]],'.',markersize=1e-10)
        ax2.set_xlim(xlim)
        ax2.set_ylim([y2[0],y2[-1]])
        labels = pylab.getp(ax2, 'yticklabels')
        pylab.setp(labels, fontsize=textsize)
        ax2.set_ylabel(yl2, fontsize=labsize)

    for jj in range(ii+1,len(ax)):
        ax[jj].set_visible(False)

    #pylab.show()

    return figg,ax


def pcolor_plot_mot1(x,y,data,time,AZ,EL,cax,xlim,ylim,xl,yl,title,text,log=0,ax=[],figg=-543,xxx=0.0,yyy=0.0,sc=1.0,textsize=8,labsize=12,xtime=1,cmap='jet'):

    pylab.ioff()

    figBG   = 'w'        # the figure background color
    axesBG  = '#f6f6f6'  # the axies background color
    figsz = (14,10)
    ncols=1; nrows=2; dx=0.1; dy=0.05
    POS=[dx,0.5,1.0-dx*2,1.0/(nrows)-dy*3]
    figg=pylab.figure(figsize=figsz, facecolor=figBG)
    ax=[]
    for aa in range(nrows):
        for bb in range(ncols):
            rect=[POS[0]+(POS[2]+dx)*bb,POS[1]-(POS[3]+dy)*aa,POS[2],POS[3]]
            ax.append(pylab.axes(rect))

    if xtime==1:
        dx=(x[-1,0]-x[0,0])/3600.0
        dx2=dx/12.0
        if dx2>0.5:
            interval=scipy.ceil(dx/12.0)
            locator = matplotlib.dates.HourLocator(interval=interval)
            formatter = matplotlib.dates.DateFormatter("%H")
        elif dx2<0.5:
            interval=scipy.ceil(dx*60.0/5.0/sc)
            locator = matplotlib.dates.MinuteLocator(interval=interval)
            formatter = matplotlib.dates.DateFormatter("%H:%M")

        x=matplotlib.dates.epoch2num(x)
        if len(xlim)!=0:
            xlim=[matplotlib.dates.epoch2num(xlim[0]),matplotlib.dates.epoch2num(xlim[-1])]
        else:
            xlim=[x[0],x[-1]]

    for ii in range(1):
        ax[ii].clear()
        if (data.shape[0]*data.shape[1])>1.0e6:
            pc=ax[ii].imshow(scipy.transpose(data),vmin=cax[0],vmax=cax[1],origin='lower',extent=(x.min(), x.max(), y.min(), y.max()),aspect='auto')
        else:
            pc=ax[ii].pcolor(scipy.transpose(x),scipy.transpose(y),scipy.transpose(data),vmin=cax[0],vmax=cax[1],cmap=pylab.get_cmap(cmap))
        if xtime==1:
            ax[ii].xaxis.set_major_locator(locator)
            ax[ii].xaxis.set_major_formatter(formatter)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        if scipy.mod(ii,ncols)==0:
            ax[ii].set_ylabel(yl, fontsize=labsize)
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        labels = pylab.getp(ax[ii], 'yticklabels')
        pylab.setp(labels, fontsize=textsize)
        ax[ii].set_xlabel(xl, fontsize=labsize)
        ax[ii].set_title(title,fontsize=labsize, horizontalalignment='center')

    pylab.axes(ax[ii])
    if log:
        cl=pylab.colorbar(pc,orientation='vertical',format=pylab.FormatStrFormatter('$10^{%.1f}$'))
    else:
        cl=pylab.colorbar(pc,orientation='vertical')
    cl.set_label(text,fontsize=labsize*1.25)

    try:
        tmp=ax[ii].get_position()
        gp=scipy.array([tmp.xmin,tmp.ymin,tmp.width,tmp.height])
    except:
        gp=ax[ii].get_position()

    ii=ii+1
    pylab.axes(ax[ii])
    try:
        tmp=ax[ii].get_position()
        gp2=scipy.array([tmp.xmin,tmp.ymin,tmp.width,tmp.height])
    except:
        gp2=ax[ii].get_position()
    gp2[2]=gp[2]
    gp2[3]=gp[3]
    ax[ii].set_position(gp2)
    for i in range(time.shape[0]):
        t1=matplotlib.dates.epoch2num(time[i,0])
        t2=matplotlib.dates.epoch2num(time[i,1])
        ax[ii].plot([t1,t2],[AZ[i,0],AZ[i,1]],'.k-')
        ax[ii].plot([t1,t2],[EL[i,0],EL[i,1]],'.b-')
#        ax[ii].set_ylim([0.0 360.0])
        ax[ii].set_xlim(xlim)
    if xtime==1:
        ax[ii].xaxis.set_major_locator(locator)
        ax[ii].xaxis.set_major_formatter(formatter)
    labels = pylab.getp(ax[ii], 'xticklabels')
    pylab.setp(labels, fontsize=textsize)
    labels = pylab.getp(ax[ii], 'yticklabels')
    pylab.setp(labels, fontsize=textsize)
    ax[ii].set_xlabel(xl, fontsize=labsize)
    ax[ii].set_ylabel('Degrees', fontsize=labsize)
    ax[ii].legend( ('Azimuth', 'Elevation') )

    for jj in range(ii+1,len(ax)):
        ax[jj].set_visible(False)

    return figg,ax

def pcolor_plot(x,y,data,cax,xlim,ylim,xl,yl,title,text,bmcodes,log=0,ax=[],figg=-543,xxx=0.0,yyy=0.0,sc=1.0,im=0,textsize=8,labsize=12,show=1,Ibeam=[],horz=2,xtime=1,doResCol=0,maxbeams=30,cmap='jet',ncols=-1,nrows=-1):

    pylab.ioff()

    Nbeams=data.shape[1]

    if Nbeams>maxbeams:
        Ibeam=range(maxbeams)
        Nbeams=maxbeams

    if (nrows==-1) and (ncols==-1):
        Nplots=Nbeams+1
        if Nbeams==1:
            nrows=1
            ncols=2
        else:
            nrows=int(scipy.ceil(scipy.sqrt(Nplots)))
            ncols=nrows
            if nrows*(nrows-1)>=Nplots:
                ncols=nrows-1

    if xtime==1:
        dx=(x[-1]-x[0])/3600.0
        dx2=dx/12.0
        if dx2>0.5:
            interval=scipy.ceil(dx/12.0)
            locator = matplotlib.dates.HourLocator(interval=interval)
            formatter = matplotlib.dates.DateFormatter("%H")
        elif dx2<0.5:
            interval=scipy.ceil(dx*60.0/5.0/sc)
            locator = matplotlib.dates.MinuteLocator(interval=interval)
            formatter = matplotlib.dates.DateFormatter("%H:%M")

        x=matplotlib.dates.epoch2num(x)
        if len(xlim)!=0:
            xlim=[matplotlib.dates.epoch2num(xlim[0]),matplotlib.dates.epoch2num(xlim[-1])]
        else:
            xlim=[x[0],x[-1]]

    ax0=ax
    if len(ax)==0 or pylab.gcf()!=figg:
        (figg,ax)=multi_axes(nrows,ncols)

    if len(Ibeam)==0:
        Ibeam=range(Nbeams)

    for ii in range(Nbeams):
        iiB=Ibeam[ii]
        ax[ii].clear()
        if (data.shape[0]*data.shape[2])>1.0e4:
            pc=ax[ii].imshow(scipy.transpose(data[:,iiB,:]),vmin=cax[0],vmax=cax[1],origin='lower',extent=(x.min(), x.max(), scipy.nanmin(y[iiB,:]), scipy.nanmax(y[iiB,:])),aspect='auto')
        else:
            pc=ax[ii].pcolor(x,scipy.squeeze(y[iiB,:]),scipy.transpose(data[:,iiB,:]),vmin=cax[0],vmax=cax[1],cmap=pylab.get_cmap(cmap))
        if xtime==1:
            ax[ii].xaxis.set_major_locator(locator)
            ax[ii].xaxis.set_major_formatter(formatter)
        ax[ii].set_xlim(xlim)
        ax[ii].set_ylim(ylim)
        if scipy.mod(ii,ncols)==0:
            ax[ii].set_ylabel(yl, fontsize=labsize)
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        labels = pylab.getp(ax[ii], 'yticklabels')
        pylab.setp(labels, fontsize=textsize)
        if ii>=(Nbeams-(ncols*nrows-Nbeams)):
#       if scipy.mod(ii,ncols)==0:
            ax[ii].set_xlabel(xl, fontsize=labsize)
        tt=r"$%d \ (%.1f^o \ \rm{az} , \ %.1f^o \ \rm{el})$" % (ii+1,bmcodes[iiB,1],bmcodes[iiB,2])
        ax[ii].set_title(tt, fontsize=labsize, horizontalalignment='center')
        if ii==0:
            ax[ii].text(xlim[0]+(xlim[1]-xlim[0])/2+xxx,(ylim[1]-ylim[0])*0.05*nrows+ylim[1]+yyy,title,fontsize=labsize, horizontalalignment='left')

    if horz>0:
        ii=ii+1
        if (len(ax0)==0) or (doResCol==1):
            try:
                tmp=ax[ii].get_position()
                gp=scipy.array([tmp.xmin,tmp.ymin,tmp.width,tmp.height])
            except:
                gp=ax[ii].get_position()
        #   gp[1]=gp[1]*2.
            gp[1]=gp[1]+gp[3]
            gp[3]=gp[3]/6.
            gp[1]=gp[1]-gp[3]
            ax[ii].set_position(gp)
        if log:
            cl=pylab.colorbar(pc,ax[ii],orientation='horizontal',format=pylab.FormatStrFormatter('$10^{%.1f}$'))
        else:
            cl=pylab.colorbar(pc,ax[ii],orientation='horizontal')
        ax[ii].yaxis.set_ticklabels([])
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        if horz>1:
            pylab.setp(labels,rotation='vertical')
    #   labels=scipy.linspace(cax[0],cax[1],len(labels)).tolist()
        cl.set_label(text,fontsize=labsize*1.25)
    else:
        if log:
            cl=figg.colorbar(pc,orientation='vertical',format=pylab.FormatStrFormatter('$10^{%.1f}$'))
        else:
            cl=figg.colorbar(pc,orientation='vertical')
        cl.set_label(text,fontsize=labsize*1.25)

    for jj in range(ii+1,len(ax)):
        ax[jj].set_visible(False)

    return figg,ax

def timegaps(time,data,rngOpt=[]):

    if len(rngOpt)>0:
        doRng=1
        rng2=[]

    time2=[]
    if scipy.ndim(data)==3:
        concnan=scipy.zeros((1,data.shape[1],data.shape[2]),dtype=data.dtype)*scipy.nan
    elif scipy.ndim(data)==2:
        concnan=scipy.zeros((1,data.shape[1]),dtype=data.dtype)*scipy.nan
    data2=data.copy()
    dt=scipy.median(scipy.diff(scipy.mean(time,axis=1)))
    for aa in range(time.shape[0]-1):
        time2.append(time[aa,0])
        if ( (time[aa+1,1]-time[aa,0]) > (dt*2.0) ):
            time2.append(time[aa,1])
            #print datetime.datetime.utcfromtimestamp(time[aa,1])
            if scipy.ndim(data)==3:
                data2=scipy.concatenate((data2[0:len(time2)-1,:,:],concnan,data2[len(time2)-1:,:,:]),axis=0)
            elif scipy.ndim(data)==2:
                data2=scipy.concatenate((data2[0:len(time2)-1,:],concnan,data2[len(time2)-1:,:]),axis=0)

    time2.append(time[-1,0])
    time2.append(time[-1,1])

    return scipy.array(time2), data2

def timegaps_wrng(time,data,rng):

    rng2=[]
    time2=[]
    concnan=scipy.zeros((1,data.shape[1]),dtype=data.dtype)*scipy.nan
    data2=data.copy()

    for aa in range(0,time.shape[0]-1):
        time2.append(time[aa,0])
        time2.append(time[aa,1])
        rng2.append(rng[aa,:].tolist())
        rng2.append(rng[aa,:].tolist())
        data2=scipy.concatenate((data2[0:len(time2)-1,:],concnan,data2[len(time2)-1:,:]),axis=0)

    time2.append(time[-1,0])
    time2.append(time[-1,1])
    rng2.append(rng[-1,:].tolist())
    rng2.append(rng[-1,:].tolist())

    return scipy.array(time2), data2, scipy.array(rng2)

def spc_pcolor_plot(meas,mod,fsc,alt,bmcodes,title='',figg1=-542,ax1=[],figg2=-543,ax2=[],figg3=-544,ax3=[],ylim=[100.0,500.0],clim=[0.0,1.0]):

    # compute spectra
    tmp=scipy.concatenate((meas,scipy.conjugate(meas[:,:0:-1,:])),axis=1) # hermitian extension
    smeas=scipy.real(scipy.fftpack.fftshift(scipy.fftpack.fft(tmp,axis=1),axes=[1])) # compute spectra
    smeas=scipy.swapaxes(smeas,0,1)

    # compute spectra
    tmp=scipy.concatenate((mod,scipy.conjugate(mod[:,:0:-1,:])),axis=1) # hermitian extension
    smod=scipy.real(scipy.fftpack.fftshift(scipy.fftpack.fft(tmp,axis=1),axes=[1])) # compute spectra
    smod=scipy.swapaxes(smod,0,1)

    freqs=scipy.arange(-meas.shape[1]+1,meas.shape[1])*fsc/1000.0

    DRC=0
    if figg1 != -542:
        DRC=1

    txt=r'$\rm{Power \ Spectral \ Density \ - \ Measured}$'
    I=scipy.where(smeas<0)
    dat=scipy.real(scipy.log10(smeas)); dat[I]=scipy.nan
    dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
    (figg1,ax1)=pcolor_plot(freqs,alt,dat,clim,[freqs[0],freqs[-1]],ylim,'Frequency (kHz)','Altitude (km)',title,txt,bmcodes,log=0,xtime=0,figg=figg1,ax=ax1,doResCol=DRC)

    txt=r'$\rm{Power \ Spectral \ Density \ - \ Modeled}$'
    I=scipy.where(smod<0)
    dat=scipy.real(scipy.log10(smod)); dat[I]=scipy.nan
    dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
    (figg2,ax2)=pcolor_plot(freqs,alt,dat,clim,[freqs[0],freqs[-1]],ylim,'Frequency (kHz)','Altitude (km)',title,txt,bmcodes,log=0,xtime=0,figg=figg2,ax=ax2,doResCol=DRC)

    txt=r'$\rm{Power \ Spectral \ Density \ - \ Residual}$'
    I=scipy.where(smod<1e-30)
    dat=(smod-smeas)/scipy.repeat(scipy.sum(smod,axis=0)[scipy.newaxis],smod.shape[0],axis=0); dat[I]=scipy.nan
    dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
    (figg3,ax3)=pcolor_plot(freqs,alt,dat,[-0.02,0.02],[freqs[0],freqs[-1]],ylim,'Frequency (kHz)','Altitude (km)',title,txt,bmcodes,log=0,xtime=0,figg=figg3,ax=ax3,doResCol=DRC)

    return figg1,ax1,figg2,ax2,figg3,ax3

def pcolor_plot_all(RF,clims=[[10,12],[0,1500],[0,3000],[0,4],[-500,500]],ylim=[],tlim=[],txMax=1.0e6,mi=[]):
    if len(ylim)==0:
        ylim=[RF.NePower['Altitude'][scipy.where(scipy.isfinite(RF.NePower['Altitude']))].min()/1000,RF.NePower['Altitude'][scipy.where(scipy.isfinite(RF.NePower['Altitude']))].max()/1000]
    if RF.FITOPTS['DO_FITS']:
        ylim0=[]
        if len(ylim0)==0:
            ylim0=[RF.FITS['Altitude'][scipy.where(scipy.isfinite(RF.FITS['Altitude']))].min()/1000,RF.FITS['Altitude'][scipy.where(scipy.isfinite(RF.FITS['Altitude']))].max()/1000]

    ttime=RF.Time['UnixTime']
    title= "%d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (RF.Time['Month'][0,0],RF.Time['Day'][0,0],RF.Time['Year'][0,0],RF.Time['dtime'][0,0],RF.Time['Month'][-1,1],RF.Time['Day'][-1,1],RF.Time['Year'][-1,1],RF.Time['dtime'][-1,1])

    Ttotal=(ttime[-1,-1]-ttime[0,0])/3600.0
    Ntrecs=scipy.ceil(Ttotal/txMax)
    print Ntrecs

    iStart=0;
    for iTime in range(Ntrecs):
        iEnd=scipy.where(ttime[:,-1]<=(ttime[iStart,0]+txMax*3600.0))[0]
        iEnd=iEnd[-1]
        tlim=[iStart,iEnd]
        iStart=iEnd+1
        if Ntrecs>1:
            txtra=' ' + str(iTime)
        else:
            txtra=''

        title2= "%d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (RF.Time['Month'][tlim[0],0],RF.Time['Day'][tlim[0],0],RF.Time['Year'][tlim[0],0],RF.Time['dtime'][tlim[0],0],RF.Time['Month'][tlim[-1],1],RF.Time['Day'][tlim[-1],1],RF.Time['Year'][tlim[-1],1],RF.Time['dtime'][tlim[-1],1])

        dt=RF.NePower['Ne_NoTr'][tlim[0]:(tlim[1]+1)]
        I=scipy.where(dt<0)

        rng=scipy.squeeze(RF.NePower['Range'][0,:])
        ytxt='Range (km)'

        xlim=[ttime[tlim[0],0],ttime[tlim[1],1]]
        tttime=ttime[tlim[0]:(tlim[1]+1)]

        # plot density, No Te/Ti
        clim=clims[0]
        txt=r'$\rm{Ne - no Tr} \ (\rm{m}^{-3})$'
    #   dat=scipy.real(RF.NePower['Ne_NoTr'])
        dat=scipy.real(scipy.log10(RF.NePower['Ne_NoTr']))[tlim[0]:(tlim[1]+1)]
        dat[I]=scipy.nan
        time,dat=timegaps(tttime,dat)
        dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
        (figg,ax)=pcolor_plot(time,RF.NePower['Altitude']/1000,dat,clim,xlim,ylim,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES,log=1)
        if RF.OPTS['saveplots']==1:
            if os.path.exists(RF.OPTS['plotsdir']):
                oname=title + '_NePower_NoTr' + txtra + '.png'
                figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

        # plot density, SNR
        clim=[-20.0,10.0]
        txt=r'$\rm{SNR} \ (\rm{dB})$'
        dat=10.0*scipy.real(scipy.log10(RF.NePower['SNR']))[tlim[0]:(tlim[1]+1)]
        dat[I]=scipy.nan
        time,dat=timegaps(tttime,dat)
        dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
        (figg,ax)=pcolor_plot(time,RF.NePower['Altitude']/1000,dat,clim,xlim,ylim,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES,log=0)
        if RF.OPTS['saveplots']==1:
            if os.path.exists(RF.OPTS['plotsdir']):
                oname=title + '_SNR' + txtra + '.png'
                figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

        if RF.FITOPTS['DO_FITS']:

            dt=RF.FITS['Ne'][tlim[0]:(tlim[1]+1)]
            I=scipy.where(dt<0)

            NION = RF.FITS['Fits'].shape[3]-1
            if len(mi)==0:
                mi=scipy.zeros((NION))

            # plot density
            clim=clims[0]
            txt=r'$\rm{Ne} \ (\rm{m}^{-3})$'
            dat=scipy.real(scipy.log10(RF.FITS['Ne']))[tlim[0]:(tlim[1]+1)]
            dat[I]=scipy.nan
            time,dat=timegaps(tttime,dat)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            print time.shape
            print RF.FITS['Altitude'].shape
            print dat.shape
            (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES,log=1)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Ne' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot ion temp
            clim=clims[1]
            for ion in range(NION):
                txt=r'$\rm{Ti} \ (\rm{K}) \ - %d \ {\rm amu}$' % (mi[ion])
                dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,1]
                dat[I]=scipy.nan
                time,dat=timegaps(tttime,dat)
                dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
                (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES)
                if RF.OPTS['saveplots']==1:
                    if os.path.exists(RF.OPTS['plotsdir']):
                        oname=title + '_Ti'+txtra + '-' + str(ion) + '.png'
                        figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot elec temp
            clim=clims[2]
            txt=r'$\rm{Te} \ (\rm{K})$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,-1,1]
            dat[I]=scipy.nan
            time,dat=timegaps(tttime,dat)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Te' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot te/ti
            clim=clims[3]
            for ion in range(NION):
                txt=r'$\rm{Te/Ti} \ - %d \ {\rm amu}$' % (mi[ion])
                dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,-1,1]/RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,1]
                dat[I]=scipy.nan
                time,dat=timegaps(tttime,dat)
                dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
                (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES)
                if RF.OPTS['saveplots']==1:
                    if os.path.exists(RF.OPTS['plotsdir']):
                        oname=title + '_Tr' +txtra + '-' + str(ion) + '.png'
                        figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot M+
            if len(clims)>5:
                clim=clims[5]
            else:
                clim=[0,1]
            for ion in range(NION):
                txt=r'$\rm{Ion \ fraction} \ - %d \ {\rm amu}$' % (mi[ion])
                dat=scipy.ma.masked_where(dt<0,RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,0])
                dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,0]
                dat[I]=scipy.nan
                time,dat=timegaps(tttime,dat)
                dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
                (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES)
                if RF.OPTS['saveplots']==1:
                    if os.path.exists(RF.OPTS['plotsdir']):
                        oname=title + '_IonFrac' +txtra + '-' + str(ion) + '.png'
                        figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot nu_in
            clim=[-2,5]
            for ion in range(NION):
                txt=r'$\rm{\nu_{in}} \ \ (\rm{s^{-1}}) \ - %d \ {\rm amu}$' % (mi[ion])
                dat=scipy.log10(RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,2])
                dat[I]=scipy.nan
                time,dat=timegaps(tttime,dat)
                dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
                (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES,log=1)
                if RF.OPTS['saveplots']==1:
                    if os.path.exists(RF.OPTS['plotsdir']):
                        oname=title + '_nuin' + txtra + '-' + str(ion) + '.png'
                        figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))


            # plot drift
            clim=clims[4]
            for ion in range(NION):
                txt=r'$\rm{Vlos} \ (\rm{m/s}) \ - %d \ {\rm amu}$' % (mi[ion])
                dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,ion,3]
                dat[I]=scipy.nan
                time,dat=timegaps(tttime,dat)
                dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
                (figg,ax)=pcolor_plot(time,RF.FITS['Altitude']/1000,dat,clim,xlim,ylim0,'Time (UT)','Altitude (km)',title2,txt,RF.BMCODES)
                if RF.OPTS['saveplots']==1:
                    if os.path.exists(RF.OPTS['plotsdir']):
                        oname=title + '_Vlos' + txtra + '-' + str(ion) + '.png'
                        figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            pylab.close('all')

    return figg,ax

def pcolor_plot_all_mot1(RF,clims=[[10,12],[0,1500],[0,3000],[0,4],[-500,500]],ylim=[],tlim=[],ylim0=[],doAlt=1,txMax=24.0,mi=[]):

    if len(ylim)==0:
        ylim=[RF.NePower['Range'][scipy.where(scipy.isfinite(RF.NePower['Range']))].min()/1000,RF.NePower['Range'][scipy.where(scipy.isfinite(RF.NePower['Range']))].max()/1000]

    ttime=RF.Time['UnixTime']

    title= "%d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (RF.Time['Month'][0,0],RF.Time['Day'][0,0],RF.Time['Year'][0,0],RF.Time['dtime'][0,0],RF.Time['Month'][-1,1],RF.Time['Day'][-1,1],RF.Time['Year'][-1,1],RF.Time['dtime'][-1,1])

    Ttotal=(ttime[-1,-1]-ttime[0,0])/3600.0
    Ntrecs=scipy.ceil(Ttotal/txMax)

    iStart=0;
    for iTime in range(Ntrecs):
        iEnd=scipy.where(ttime[:,-1]<=(ttime[iStart,0]+txMax*3600.0))[0]
        iEnd=iEnd[-1]
        tlim=[iStart,iEnd]
        iStart=iEnd+1
        if Ntrecs>1:
            txtra=' ' + str(iTime)
        else:
            txtra=''

        title2= "%d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (RF.Time['Month'][tlim[0],0],RF.Time['Day'][tlim[0],0],RF.Time['Year'][tlim[0],0],RF.Time['dtime'][tlim[0],0],RF.Time['Month'][tlim[-1],1],RF.Time['Day'][tlim[-1],1],RF.Time['Year'][tlim[-1],1],RF.Time['dtime'][tlim[-1],1])

        dt=RF.NePower['Ne_NoTr'][tlim[0]:(tlim[1]+1)]
        I=scipy.where(dt<0)

        AZ=RF.Antenna['Azimuth'][tlim[0]:(tlim[1]+1)]
        EL=RF.Antenna['Elevation'][tlim[0]:(tlim[1]+1)]

        rng=scipy.squeeze(RF.NePower['Range'][0,:])
        ytxt='Range (km)'

        xlim=[ttime[tlim[0],0],ttime[tlim[1],1]]
        tttime=ttime[tlim[0]:(tlim[1]+1)]

        # plot density, No Te/Ti
        clim=clims[0]
        txt=r'$\rm{Ne - no Tr} \ (\rm{m}^{-3})$'
        dat=scipy.real(scipy.log10(RF.NePower['Ne_NoTr']))[tlim[0]:(tlim[1]+1)]
        dat[I]=scipy.nan
        time,dat=timegaps(tttime,scipy.squeeze(dat))
        dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
        (figg,ax)=pcolor_plot_mot1(time[:,scipy.newaxis],rng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=1)
        if RF.OPTS['saveplots']==1:
            if os.path.exists(RF.OPTS['plotsdir']):
                oname=title + '_NePower_NoTr' + txtra + '.png'
                figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

        # plot density, SNR
        clim=[-20.0,10.0]
        txt=r'$\rm{SNR} \ (\rm{dB})$'
        dat=10.0*scipy.real(scipy.log10(RF.NePower['SNR']))[tlim[0]:(tlim[1]+1)]
        dat[I]=scipy.nan
        time,dat=timegaps(tttime,scipy.squeeze(dat))
        dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
        (figg,ax)=pcolor_plot_mot1(time[:,scipy.newaxis],rng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
        if RF.OPTS['saveplots']==1:
            if os.path.exists(RF.OPTS['plotsdir']):
                oname=title + '_NePower_SNR' + txtra + '.png'
                figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

        pylab.close('all')

        if RF.FITOPTS['DO_FITS']:

            alt=scipy.squeeze(RF.FITS['Altitude'][tlim[0]:(tlim[1]+1)])
            rng=scipy.squeeze(RF.FITS['Range'][tlim[0]:(tlim[1]+1)])

            dt=RF.FITS['Ne'][tlim[0]:(tlim[1]+1)]
            I=scipy.where(dt<0)

            if doAlt==1:
                dr=scipy.diff(alt,axis=1)
                dr=scipy.concatenate((dr,dr[:,[-1]]),axis=1)/2.0
                rng=scipy.concatenate((alt-dr,alt[:,[-1]]+dr[:,[-1]]),axis=1)
                ytxt='Altitude (km)'
                if len(ylim0)==0:
                    ylim=[alt[scipy.where(scipy.isfinite(alt))].min()/1000,alt[scipy.where(scipy.isfinite(alt))].max()/1000]
            else:
                dr=scipy.diff(rng,axis=1)
                dr=scipy.concatenate((dr,dr[:,[-1]]),axis=1)/2.0
                rng=scipy.concatenate((rng-dr,rng[:,[-1]]+dr[:,[-1]]),axis=1)
                ytxt='Range (km)'
                if len(ylim0)==0:
                    ylim=[rng[scipy.where(scipy.isfinite(rng))].min()/1000,rng[scipy.where(scipy.isfinite(rng))].max()/1000]

            # plot density
            clim=clims[0]
            txt=r'$\rm{Ne} \ (\rm{m}^{-3})$'
            dat=scipy.real(scipy.log10(RF.FITS['Ne']))[tlim[0]:(tlim[1]+1)]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=1)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Ne' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot ion temp
            clim=clims[1]
            txt=r'$\rm{Ti} \ (\rm{K})$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,0,1]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Ti' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot elec temp
            clim=clims[2]
            txt=r'$\rm{Te} \ (\rm{K})$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,-1,1]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Te' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot te/ti
            clim=clims[3]
            txt=r'$\rm{Te/Ti}$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,-1,1]/RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,0,1]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Tr' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot drift
            clim=clims[4]
            txt=r'$\rm{Vlos} \ (\rm{m/s})$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,0,3]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_Vlos' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot M+
            if len(clims)>5:
                clim=clims[5]
            else:
                clim=[0,1]
            txt=r'$\rm{O^+ frac}$'
            dat=RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,0,0]
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=0)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_OpFrac' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

            # plot nu_in
            clim=[-2,5]
            txt=r'$\rm{\nu_{in}} \ \ (\rm{s^{-1}})$'
            dat=scipy.log10(RF.FITS['Fits'][tlim[0]:(tlim[1]+1),:,:,0,2])
            dat[I]=scipy.nan
            time,dat,trng=timegaps_wrng(tttime,scipy.squeeze(dat),rng)
            time=scipy.repeat(time[:,scipy.newaxis],dat.shape[1]+1,axis=1)
            dat=scipy.ma.masked_where(scipy.isnan(dat),dat)
            (figg,ax)=pcolor_plot_mot1(time,trng/1000,dat,tttime,AZ,EL,clim,xlim,ylim,'Time (UT)',ytxt,title2,txt,log=1)
            if RF.OPTS['saveplots']==1:
                if os.path.exists(RF.OPTS['plotsdir']):
                    oname=title + '_nuin' + txtra + '.png'
                    figg.savefig(os.path.join(RF.OPTS['plotsdir'],oname))

    return figg,ax

def replot_pcolor_all(fname,saveplots=0,opath='.',clims=[[10,12],[0,1500],[0,3000],[0,4],[-500,500]],ylim=[],tlim=[],txMax=1.0e6):

    RF=MyClass()

    h5file=tables.openFile(fname)
    output={}
    for array in h5file.listNodes('/',classname = 'Array'):
        output[array.name]=array.read()
    for group in h5file.walkGroups("/"):
        output[group._v_pathname]={}
        for array in h5file.listNodes(group, classname = 'Array'):
            output[group._v_pathname][array.name]=array.read()
    h5file.close()

    RF.Time={}
    RF.NePower={}
    RF.OPTS={}
    RF.FITS={}
    RF.VVELS={}
    RF.Antenna={}
#   RF.DO_VVELS=0
    RF.FITOPTS={}
    RF.FITOPTS['DO_FITS']=0
    mi=[]

    RF.OPTS['saveplots']=saveplots
    RF.FITTER_PATH=''
    RF.OPTS['outputpath']=''
    RF.OPTS['plotsdir']=opath

    RF.Time['UnixTime']=output['/Time']['UnixTime']
    RF.Time['dtime']=output['/Time']['dtime']
    RF.Time['Month']=output['/Time']['Month']
    RF.Time['Day']=output['/Time']['Day']
    RF.Time['Year']=output['/Time']['Year']

    if output.has_key('/Antenna'):
        RF.MOTION_TYPE=1
        RF.Antenna['Azimuth']=output['/Antenna']['Azimuth']
        RF.Antenna['Elevation']=output['/Antenna']['Elevation']
    else:
        RF.MOTION_TYPE=0

    if output.has_key('/FittedParams'):
        RF.FITOPTS['DO_FITS']=1
        RF.FITS['Fits']=output['/FittedParams']['Fits']
        RF.FITS['Ne']=output['/FittedParams']['Ne']
        RF.FITS['Altitude']=output['/FittedParams']['Altitude']
        RF.FITS['Range']=output['/FittedParams']['Range']
        mi=output['/FittedParams']['IonMass']

    RF.NePower['Ne_NoTr']=output['/NeFromPower']['Ne_NoTr']
    RF.NePower['Ne_Mod']=output['/NeFromPower']['Ne_Mod']
    RF.NePower['Altitude']=output['/NeFromPower']['Altitude']
    RF.NePower['Range']=output['/NeFromPower']['Range']
    RF.NePower['SNR']=output['/NeFromPower']['SNR']

#   if output.has_key('/VectorVels'):
#       RF.DO_VVELS=1
#       RF.VVELS['Vest']=output['/VectorVels']['Vest']
#       RF.VVELS['dVest']=output['/VectorVels']['dVest']
#       RF.VVELS['plat']=output['/VectorVels']['plat']

    if RF.MOTION_TYPE==0:
        RF.BMCODES=output['BeamCodes']
        pcolor_plot_all(RF,clims=clims,ylim=ylim,tlim=tlim,txMax=txMax,mi=mi)
    elif RF.MOTION_TYPE==1:
        pcolor_plot_all_mot1(RF,clims=clims,ylim=ylim,tlim=tlim,txMax=txMax,mi=mi)

    return RF

def acf_plot(meas,errs,mod,alt,bmcodes,title='',ax=[],figg=-542,maxbeams=10,Ibeams=[],maxlag=-1):

    pylab.ioff()

    textsize = 8        # size for axes text
    labsize = 10

    (Nbeams,Nlags,Nhts)=meas.shape
    nrows=1
    if len(Ibeams)==0:
        ncols=Nbeams
        Ibeams=range(ncols)
    else:
        ncols=len(Ibeams)
    if ncols>maxbeams:
        ncols=maxbeams
        Ibeams=range(ncols)
    ii=scipy.where(scipy.isfinite(alt))
    ylim=[alt[ii].min()-25,alt[ii].max()+25]
    xl='Lag'
    yl='Altitude (km)'

    if maxlag==-1:
        maxlag=Nlags

    if len(ax)==0 or figg==-542:
        (figg,ax)=multi_axes(nrows,ncols)
    pylab.figure(figg.number)

    for ii in range(ncols):
        rr=Ibeams[ii]
        ax[ii].clear()
        ax[ii].hold(True)
        for jj in range(Nhts):
            try:
                dh=(alt[rr,jj+1]-alt[rr,jj])/2.0
            except:
                ''
            Iy=scipy.where(scipy.isfinite(mod[rr,:,jj]))[0]
            if jj==0:
                ml=len(Iy)
            elif len(Iy)>ml:
                ml=len(Iy)
            if len(Iy)>0:
                sc=1./scipy.nanmax(scipy.absolute(meas[rr,Iy,jj]))*dh
                ax[ii].errorbar(scipy.arange(len(Iy)),alt[rr,jj]+sc*meas[rr,Iy,jj].imag,yerr=sc*scipy.sqrt(errs[rr,Iy,jj]),fmt='r-')
                ax[ii].errorbar(scipy.arange(len(Iy)),alt[rr,jj]+sc*meas[rr,Iy,jj].real,yerr=sc*scipy.sqrt(errs[rr,Iy,jj]),fmt='b-')
                ax[ii].plot(range(len(Iy)),alt[rr,jj]+sc*mod[rr,Iy,jj].real,'k')
                ax[ii].plot(range(len(Iy)),alt[rr,jj]+sc*mod[rr,Iy,jj].imag,'k')
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlim([0,ml])
        ax[ii].set_xlabel(xl, fontsize=labsize)
        if ii==0:
            ax[ii].set_ylabel(yl, fontsize=labsize)
            ax[ii].text(1,(ylim[1]-ylim[0])*0.1+ylim[1],title,fontsize=labsize, horizontalalignment='left')
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        labels = pylab.getp(ax[ii], 'yticklabels')
        pylab.setp(labels, fontsize=textsize)

        tt=r"$(%.1f^o \ \rm{az} , \ %.1f^o \ \rm{el})$" % (bmcodes[rr,1],bmcodes[rr,2])
        ax[ii].set_title(tt, fontsize=labsize, horizontalalignment='center')

    #pylab.show()

    return figg,ax

def spc_plot(meas,errs,mod,alt,bmcodes,title='',ax=[],figg=-542,maxbeams=10,Ibeams=[],measSpc=0,fsc=1.0,Ialt=-1,ylim=-1,scfac=1.0,ploterrs=1,zeroline=1):

    pylab.ioff()

    textsize = 8        # size for axes text
    labsize = 10

    (Nbeams,Nlags,Nhts)=meas.shape
    nrows=1
    if len(Ibeams)==0:
        ncols=Nbeams
        Ibeams=range(ncols)
    else:
        ncols=len(Ibeams)
    if ncols>maxbeams:
        ncols=maxbeams
        Ibeams=range(ncols)
    ii=scipy.where(scipy.isfinite(alt))
    if ylim==-1:
        ylim=[alt[ii].min()-25,alt[ii].max()+25]
    xl='Frequency'
    yl='Altitude (km)'
    dh=(alt[Ibeams[0],1]-alt[Ibeams[0],0])/2.0

    if measSpc==0:
        Iy=scipy.where(scipy.absolute(scipy.nansum(mod[0,:,:],1))>0.0)[0]
        meas=meas[:,Iy,:]; mod=mod[:,Iy,:]; errs=errs[:,Iy,:]
        if scipy.all(scipy.isnan(scipy.absolute(meas[:,0,:]))):
            meas[:,0,:]=mod[:,0,:]

        # compute spectra
        tmp=scipy.concatenate((meas,scipy.conjugate(meas[:,:0:-1,:])),axis=1) # hermitian extension
        smeas=scipy.real(scipy.fftpack.fftshift(scipy.fftpack.fft(tmp,axis=1),axes=[1])) # compute spectra

        # compute spectra
        tmp=scipy.concatenate((mod,scipy.conjugate(mod[:,:0:-1,:])),axis=1) # hermitian extension
        smod=scipy.real(scipy.fftpack.fftshift(scipy.fftpack.fft(tmp,axis=1),axes=[1])) # compute spectra

        freqs=scipy.arange(-meas.shape[1]+1,meas.shape[1])*fsc

    else:
        smeas=meas; smod=mod
        freqs=scipy.arange(-scipy.ceil(meas.shape[1]/2.0)+1,scipy.ceil(meas.shape[1]/2.0))*fsc

    if len(ax)==0 or figg==-542:
        (figg,ax)=multi_axes(nrows,ncols)
    pylab.figure(figg.number)

    for ii in range(ncols):
        rr=Ibeams[ii]
        ax[ii].clear()
        ax[ii].hold(True)
        if Ialt==-1:
            Ialtt=range(Nhts)
        else:
            Ialtt=Ialt[ii]
        for jj in Ialtt:
            try:
                dh=(alt[rr,jj+1]-alt[rr,jj])/2.0
            except:
                ''
            sc=scfac/scipy.absolute(smeas[rr,:,jj]).max()*dh
            derr=scipy.sqrt(errs[rr,1,jj])/meas[rr,1,jj]
            if ploterrs:
                ax[ii].errorbar(freqs,alt[rr,jj]+sc*smeas[rr,:,jj],yerr=sc*smeas[rr,:,jj]*derr,fmt='r-')
            else:
                ax[ii].plot(freqs,alt[rr,jj]+sc*smeas[rr,:,jj],'r-')
            ax[ii].plot(freqs,alt[rr,jj]+sc*smod[rr,:,jj],'k')
            ax[ii].plot(freqs,freqs/freqs*alt[rr,jj],'k--')
        ax[ii].set_ylim(ylim)
        ax[ii].set_xlim([freqs[0],freqs[-1]])
        ax[ii].set_xlabel(xl, fontsize=labsize)
        if ii==0:
            ax[ii].set_ylabel(yl, fontsize=labsize)
            ax[ii].text(1,(ylim[1]-ylim[0])*0.1+ylim[1],title,fontsize=labsize, horizontalalignment='left')
        labels = pylab.getp(ax[ii], 'xticklabels')
        pylab.setp(labels, fontsize=textsize)
        labels = pylab.getp(ax[ii], 'yticklabels')
        pylab.setp(labels, fontsize=textsize)

        tt=r"$(%.1f^o \ \rm{az} , \ %.1f^o \ \rm{el})$" % (bmcodes[rr,1],bmcodes[rr,2])
        ax[ii].set_title(tt, fontsize=labsize, horizontalalignment='center')

    #pylab.show()

    return figg,ax

def test_plot(RF,irec,ax=[],figg=-542,xlims=[(0.01,10),(0,3),(-1,1)],maxbeams=10,Ibeams=[],dofrac=0):

    nrows=3
    if dofrac:
        nrows=nrows+1

    if len(Ibeams)==0:
        ncols=RF.Nbeams
        Ibeams=range(ncols)
    else:
        ncols=len(Ibeams)
    nalts=RF.FITS['Altitude'].shape[1]

    if ncols>maxbeams:
        ncols=maxbeams
        Ibeams=range(ncols)

    #I=scipy.where(RF.FITS['FitInfo']['fitcode']<=0)

    # densities
    in1a=scipy.reshape(RF.FITS['Ne'][Ibeams,:],(ncols,nalts,1))/1.0e11
    in1b=scipy.reshape(RF.FITS['Altitude'][Ibeams,:],(ncols,nalts,1))/1000.0
    in1c=scipy.reshape(RF.FITS['dNe'][Ibeams,:],(ncols,nalts,1))/1.0e11

    # temperatures
    in2a=RF.FITS['Fits'][Ibeams,:,:,1]/1000.0
    in2b=scipy.repeat(scipy.reshape(RF.FITS['Altitude'][Ibeams,:],(ncols,nalts,1)),RF.FITS['Fits'].shape[2],axis=2)/1000.0
    in2c=RF.FITS['Errors'][Ibeams,:,:,1]/1000.0

    # velocities
#    in3a=scipy.reshape(RF.FITS['Fits'][Ibeams,:,0,3],(ncols,nalts,1))/100.0
#    in3b=scipy.reshape(RF.FITS['Altitude'][Ibeams,:],(ncols,nalts,1))/1000.0
#    in3c=scipy.reshape(RF.FITS['Errors'][Ibeams,:,0,3],(ncols,nalts,1))/100.0
    in3a=RF.FITS['Fits'][Ibeams,:,:,3]/100.0
    in3b=scipy.repeat(scipy.reshape(RF.FITS['Altitude'][Ibeams,:],(ncols,nalts,1)),RF.FITS['Fits'].shape[2],axis=2)/1000.0
    in3c=RF.FITS['Errors'][Ibeams,:,:,3]/100.0


    # fractions
    in4a=RF.FITS['Fits'][Ibeams,:,:,0]
    in4b=scipy.repeat(scipy.reshape(RF.FITS['Altitude'][Ibeams,:],(ncols,nalts,1)),RF.FITS['Fits'].shape[2],axis=2)/1000.0
    in4c=RF.FITS['Errors'][Ibeams,:,:,0]

    title= "%d-%d-%d %.3f-%.3f UT" % (RF.Time['Month'][0],RF.Time['Day'][0],RF.Time['Year'][0],RF.Time['dtime'][0],RF.Time['dtime'][1])
    xlabels=[r'$\rm{Ne} \ (10^{11} \ \rm{m}^{-3})$', r'$\rm{T} \ (10^3 \ \rm{K})$', r'$\rm{Vlos} \ (10^2 \ \rm{m/s})$']
    ylabels=['Altitude (km)','Altitude (km)','Altitude (km)','Altitude (km)']

    if dofrac:
        xlabels.append(r'$\rm{Fraction}$')
        ylabels.append('Altitude (km)')
        if len(xlims)==3:
            xlims.append((0,1))

        if len(ax)==0 or figg==-542:
            (figg,ax)=multi_axes(nrows,ncols,dx=0.015)
        pylab.figure(figg.number)
        plot_array(ax,nrows,ncols,in1a,in1b,in1c,in2a,in2b,in2c,in3a,in3b,in3c,in4a,in4b,in4c,xlims,xlabels,ylabels,title,RF.BMCODES[Ibeams,:])

    else:
        xlabels=[r'$\rm{Ne} \ (10^{11} \ \rm{m}^{-3})$', r'$\rm{T} \ (10^3 \ \rm{K})$', r'$\rm{Vlos} \ (10^2 \ \rm{m/s})$']
        ylabels=['Altitude (km)','Altitude (km)','Altitude (km)']
        title= "%d-%d-%d %.3f-%.3f UT" % (RF.Time['Month'][0],RF.Time['Day'][0],RF.Time['Year'][0],RF.Time['dtime'][0],RF.Time['dtime'][1])

        if len(ax)==0 or figg==-542:
            (figg,ax)=multi_axes(nrows,ncols,dx=0.015)
        pylab.figure(figg.number)
        plot_array(ax,nrows,ncols,in1a,in1b,in1c,in2a,in2b,in2c,in3a,in3b,in3c,xlims,xlabels,ylabels,title,RF.BMCODES[Ibeams,:])

    return figg,ax,title
