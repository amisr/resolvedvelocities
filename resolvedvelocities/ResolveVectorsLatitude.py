#!/usr/bin/env python

"""

"""


import os
import time
import datetime
import ConfigParser

import vvels
import plotVvels
import ioutils
from ioclass import OutputFileClass

import numpy as np


class vvelsLatFile(OutputFileClass):

    def __init__(self):
        self.title = 'Vector Velocities Resolved in Latitude (geographic or geomagnetic)'
        
        self.h5Paths = {
            'Params'    :   ('/ProcessingParams','Experiment Parameters'),\
            'Site'      :   ('/Site','Site Parameters'),\
            'Time'      :   ('/Time','Time Information'),\
            'VectorVels' : ('/VectorVels','Vector Velocities'),\
        }        
                        
        self.h5Attribs = {
            '/Time/dtime' : [('TITLE','Decimal Time'),('Size','Nrecords x 2 (Start and end of integration)'),('Unit','UT Hours')],\
            '/Time/UnixTime' : [('TITLE','Unix Time'),('Size','Nrecords x 2 (Start and end of integration)'),('Unit','Seconds')],\
            '/Time/MagneticLocalTime' : [('TITLE','Magnetic Local Time'),('Size','Nrecords x 2 (Start and end of integration)'),('Unit','UT Hours')],\
            '/Site/Altitude' : [('TITLE','Altitude'),('Description','Altitude of site'),('Unit','Meters')],\
            '/Site/Code' : [('TITLE','Site Code')],\
            '/Site/Latitude' : [('TITLE','Latitude'),('Description','Latitude of site'),('Unit','Degrees North')],\
            '/Site/Longitude' : [('TITLE','Longitude'),('Description','Longitude of site'),('Unit','Degrees East')],\
            '/Site/MagneticLatitude' : [('TITLE','Magnetic Latitude'),('Description','Magnetic Latitude of site'),('Unit','Degrees North')],\
            '/Site/MagneticLongitude' : [('TITLE','Magnetic Longitude'),('Description','Magnetic Longitude of site'),('Unit','Degrees East')],\
            '/Site/MagneticLocalTimeMidnight' : [('TITLE','Magnetic Local Time Midnight'),('Unit','UT Hours')],\
            '/Site/Name' : [('TITLE','Name'),('Description','Site Name')],\
            '/VectorVels/Latitude' : [('TITLE','Latitude'),('Description','Latitude of measurement'),('Unit','Degrees')],\
            '/VectorVels/MagneticLatitude' : [('TITLE','MagneticLatitude'),('Description','Geomagnetic latitude of measurement'),('Unit','Degrees')],\
            '/VectorVels/Eavg' : [('TITLE','Electric Field'),('Description','Averaged Electric Field'),('Unit','V/m')],\
            '/VectorVels/Bavg' : [('TITLE','Magnetic Field'),('Description','Averaged Magnetic Field Used for E-field Calculation'),('Unit','T')],\
            '/VectorVels/Nmeas' : [('TITLE','Number of Measurements'),('Description','Number of Measurements in Estimate')],\
            '/VectorVels/Vdir' : [('TITLE','Ion Velocity Direction'),('Description','Ion Velocity Direction'),('Unit','Degrees North of East')],\
            '/VectorVels/Vest' : [('TITLE','Ion Velocity Vector'),('Description','Ion Velocity Vector'),('Size','Nrecords x Naltitudes x 3 (Perp-North, Perp-East, Anti-Parallel)'),('Unit','m/s')],\
            '/VectorVels/Vmag' : [('TITLE','Ion Velocity Magnitude'),('Description','Ion Velocity Magnitude'),('Unit','m/s')],\
            '/VectorVels/errVdir' : [('TITLE','Error on Ion Velocity Direction'),('Description','Error on Ion Velocity Direction'),('Unit','Degrees')],\
            '/VectorVels/errVest' : [('TITLE','Error on Ion Velocity Vector'),('Description','Error on Ion Velocity Vector'),('Size','Nrecords x Naltitudes x 3 (Perp-North, Perp-East, Anti-Parallel)'),('Unit','m/s')],\
            '/VectorVels/varVest' : [('TITLE','Covariance Matrix for Ion Velocity Vector Estimate'),('Description','Ion Velocity Vector Covariance Matrix'),('Size','Nrecords x Naltitudes x 3 x 3 (Perp-North, Perp-East, Anti-Parallel)'),('Unit','(m/s)^2')],\
            '/VectorVels/errVmag' : [('TITLE','Error on Ion Velocity Magnitude'),('Description','Error on Ion Velocity Magnitude'),('Unit','m/s')],\
            '/VectorVels/errEest' : [('TITLE','Error on Electric Field'),('Description','Error on Electric Field'),('Unit','V/m')],\
            '/VectorVels/varEest' : [('TITLE','Covariance Matrix for Electric Field Estimate'),('Description','Electric Field Vector Covariance Matrix'),('Size','Nrecords x Naltitudes x 3 x 3 (Perp-North, Perp-East, Anti-Parallel)'),('Unit','(V/m)^2')],\
            '/ProcessingParams/ProcessingTimeStamp' : [('TITLE','Processing Time Stamp')],\
            '/ProcessingParams/BaudLength' : [('TITLE','Baud Length'),('Unit','Seconds')],\
            '/ProcessingParams/PulseLength' : [('TITLE','Pulse Length'),('Unit','Seconds')],\
            '/ProcessingParams/RxFrequency' : [('TITLE','Rx Frequency'),('Description','Receive frequency'),('Unit','Hertz')],\
            '/ProcessingParams/TxFrequency' : [('TITLE','Tx Frequency'),('Description','Transmit frequency'),('Unit','Hertz')],\
            '/ProcessingParams/MinAlt': [('TITLE','Min Altitude'),('Description','Minimum Altitude for E-field Calculation'),('Unit','Meters')],\
            '/ProcessingParams/MaxAlt': [('TITLE','Max Altitude'),('Description','Maximum Altitude for E-field Calculation'),('Unit','Meters')],\
            '/ProcessingParams/Covar': [('TITLE','Covariance Matrix'),('Description','A-priori Covariance Matrix (perp-North, perp-East, anti-parallel)'),('Unit','(m/s)^2')],\
            '/ProcessingParams/ErrorElim': [('TITLE','Error Elimination'),('Description','Error Elimination Criteria')],\
            '/ProcessingParams/SourceFile': [('TITLE','Source File'),('Description','Source File of Measurements')],
            '/ProcessingParams/CorrectVap': [('TITLE','Correction using Field-Aligned Velocity'),('Description','Flag indicating whether velocities were corrected using field-aligned measurements')],
            '/ProcessingParams/GeographicBinning': [('TITLE','Geographic or Geomagnetic Binning'),('Description','0: Local geomagnetic by geomagnetic latitude, 1: Geographic by geomagnetic latitude, 2: Geographic by geographic latitude')],
            '/ProcessingParams/IntegrationTime': [('TITLE','Integration Time'),('Unit','Seconds')],
            '/ProcessingParams/MinimumElectronDensity': [('TITLE','Electron Density Threshold'),('Description', 'Data below threshold is ignored'), ('Unit','e/m^-3')],
            '/ProcessingParams/VelocityOffsetCorrection': [('TITLE','Offset to LOS Velocity (Chirp)'),('Description', 'Additive velocity offset'), ('Unit','m/s')],
        }

        return
        
class vvelsLat:
    
    def __init__(self,inifiles,sec):
    
        # self.DefaultIni = '/home/asreimer/projects/AMISR/development/source_code/fitting_software/amisr-src/src/ISfit/AMISR_fitter_py/config/config_vvelsLat-default.ini'

        print inifiles
        print sec

        # parse ini file
        self.ini_parse(inifiles,sec)
        
        # input file
        self.h5file=ioutils.h5file(self.filename) # h5file instance

        # output directory
        if not os.path.exists(os.path.dirname(self.outputFname)):
            os.makedirs(os.path.dirname(self.outputFname))

        # create output names
        self.setOutPutNames()

        return
        
    def ini_parse(self,inifile,sec):        

        # setup ConfigParser object
        config=ConfigParser.ConfigParser()
        
        # # read default ini file
        # fn = config.read(self.DefaultIni)
        # if len(fn)!=1 or fn[0]!=self.DefaultIni:
        #     raise IOError, 'Unable to read default ini file %s' % self.DefaultIni

        # read supplied ini file
        config.read(inifile.split(','))

        # make sure necessary section exists
        if not config.has_section(sec):
            raise IOError, 'Configuration files must contain section %s' % sec

        self.config=config

        # set processing parameters
        # self.minAlt=eval(ioutils.ini_tool(config,sec,'minalt',required=1,defaultParm=config.get('DefaultParameters','minalt')))
        # self.maxAlt=eval(ioutils.ini_tool(config,sec,'maxalt',required=1,defaultParm=config.get('DefaultParameters','maxalt')))
        # self.ppp=eval(ioutils.ini_tool(config,sec,'ppp',required=1,defaultParm=config.get('DefaultParameters','ppp')))
        # self.covar=eval(ioutils.ini_tool(config,sec,'covar',required=1,defaultParm=config.get('DefaultParameters','covar')))
        # self.zoomWhole=eval(ioutils.ini_tool(config,sec,'zoomwhole',required=1,defaultParm=config.get('DefaultParameters','zoomwhole')))
        # self.chirp=eval(ioutils.ini_tool(config,sec,'chirp',required=1,defaultParm=config.get('DefaultParameters','chirp')))
        # self.neMin=eval(ioutils.ini_tool(config,sec,'nemin',required=1,defaultParm=config.get('DefaultParameters','nemin')))
        # self.byGeo=eval(ioutils.ini_tool(config,sec,'bygeo',required=1,defaultParm=config.get('DefaultParameters','bygeo')))
        # self.Time2Integrate=eval(ioutils.ini_tool(config,sec,'time2integrate',required=1,defaultParm=config.get('DefaultParameters','time2integrate')))
        # self.plats=eval(ioutils.ini_tool(config,sec,'plats',required=1,defaultParm=config.get('DefaultParameters','plats')))
        # self.CorrectVap=eval(ioutils.ini_tool(config,sec,'correctvap',required=1,defaultParm=config.get('DefaultParameters','correctvap')))
        # self.upBcode=eval(ioutils.ini_tool(config,sec,'upbcode',required=1,defaultParm=config.get('DefaultParameters','upbcode')))
        # self.beams2use=eval(ioutils.ini_tool(config,sec,'beams2use',required=1,defaultParm=list()))
        # # set plotting parameters
        # self.makeplot=eval(ioutils.ini_tool(config,sec,'makeplot',required=1,defaultParm=config.get('DefaultParameters','makeplot')))
        # self.clim=eval(ioutils.ini_tool(config,sec,'clim',required=1,defaultParm=config.get('DefaultParameters','clim')))
        # self.sc=eval(ioutils.ini_tool(config,sec,'sc',required=1,defaultParm=config.get('DefaultParameters','sc')))

        self.minAlt=eval(ioutils.ini_tool(config,sec,'minalt',required=1))
        self.maxAlt=eval(ioutils.ini_tool(config,sec,'maxalt',required=1))
        self.ppp=eval(ioutils.ini_tool(config,sec,'ppp',required=1))
        self.covar=eval(ioutils.ini_tool(config,sec,'covar',required=1))
        self.zoomWhole=eval(ioutils.ini_tool(config,sec,'zoomwhole',required=1))
        self.chirp=eval(ioutils.ini_tool(config,sec,'chirp',required=1))
        self.neMin=eval(ioutils.ini_tool(config,sec,'nemin',required=1))
        self.byGeo=eval(ioutils.ini_tool(config,sec,'bygeo',required=1))
        self.Time2Integrate=eval(ioutils.ini_tool(config,sec,'time2integrate',required=1))
        self.plats=eval(ioutils.ini_tool(config,sec,'plats',required=1))
        self.CorrectVap=eval(ioutils.ini_tool(config,sec,'correctvap',required=1))
        self.upBcode=eval(ioutils.ini_tool(config,sec,'upbcode',required=1))
        self.beams2use=eval(ioutils.ini_tool(config,sec,'beams2use',required=1))
        # set plotting parameters
        self.makeplot=eval(ioutils.ini_tool(config,sec,'makeplot',required=1))
        self.clim=eval(ioutils.ini_tool(config,sec,'clim',required=1))
        self.sc=eval(ioutils.ini_tool(config,sec,'sc',required=1))

        
        # get filenames
        self.saveout=1
        self.filename=ioutils.ini_tool(config,sec,'inputFilename',required=1)
        self.outputFname=ioutils.ini_tool(config,sec,'outputFilename',required=1)

        return
        
    def setOutDicts(self,dat):

        self.Time = {
            'UnixTime'  :   None,\
            'dtime' :   None,\
            'MagneticLocalTime' : None
        }
        
        self.Site = dat['/Site']
        
        self.VectorVels = {
            'Latitude' : None,\
            'MagneticLatitude' : None,\
            'Nmeas' : None,\
            'Vdir' : None,\
            'Vest' : None,\
            'Vmag' : None,\
            'Edir' : None,\
            'Eest' : None,\
            'Emag' : None,\
            'errVmag' : None,\
            'errVest' : None,\
            'errVdir' : None,\
            'errEmag' : None,\
            'errEest' : None,\
            'errEdir' : None,\
        }
        
        self.Params = {
            'MinAlt'    :   self.minAlt*1e3, \
            'MaxAlt'    :   self.maxAlt*1e3, \
            'Covar'     :   self.covar, \
            'ErrorElim' :   self.ppp, \
            'SourceFile'    :   os.path.basename(self.filename), \
            'PulseLength'   :   None ,\
            'BaudLength'   :   None ,\
            'RxFrequency'   :   None ,\
            'TxFrequency'   :   None ,\
            'ProcessingTimeStamp'    : time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),\
            'TitleString'   :   self.titleString,\
            'CorrectVap'    :   self. CorrectVap,\
            'GeographicBinning' :   self.byGeo,\
            'IntegrationTime'   :   self.Time2Integrate,\
            'MinimumElectronDensity' : self.neMin,\
            'VelocityOffsetCorrection': self.chirp,\
        }

        try: self.Params['PulseLength']=dat['/ProcessingParams']['PulseLength']
        except: self.Params['PulseLength']=0.0
        try: self.Params['BaudLength']=dat['/ProcessingParams']['BaudLength']
        except: self.Params['BaudLength']=0.0
        try: self.Params['RxFrequency']=dat['/ProcessingParams']['RxFrequency']
        except: self.Params['RxFrequency']=0.0
        try: self.Params['TxFrequency']=dat['/ProcessingParams']['TxFrequency']
        except: self.Params['TxFrequency']=0.0
    
        return
        
    def createOutputFile(self,outputFname):
        
        # create output file
        ofile = vvelsLatFile()
        ofile.createFile(outputFname+'.h5')

        # create groups
        ofile.createh5groups()
        
        # create arrays
        ofile.createStaticArray(ofile.h5Paths['Params'][0],self.Params,keys2do=self.Params.keys())
        ofile.createStaticArray(ofile.h5Paths['Time'][0],self.Time,keys2do=self.Time.keys())
        ofile.createStaticArray(ofile.h5Paths['VectorVels'][0],self.VectorVels,keys2do=self.VectorVels.keys())
        ofile.createStaticArray(ofile.h5Paths['Site'][0],self.Site,keys2do=self.Site.keys())
        
        # set attributes in output file
        ofile.setAtrributes()
        
        # close file
        ofile.closeFile()
        
        return
           
    def replotFromOutput(self,sc=None,clim=None,createAll=0,binByDay=1,txMax=1.0e6):
    
        if sc != None:
            self.sc = sc
        if clim != None:
            self.clim = clim
        
        
        
        for ofname in self.outputFnames:
            ofid = ioutils.h5file(ofname+'.h5')
            dat = ofid.readWholeh5file()
            self.VectorVels = dat['/VectorVels']
            self.Time = dat['/Time']
        
            timeout = self.Time['UnixTime']
            try:
                self.titleString
            except:
                ds = datetime.datetime.utcfromtimestamp(timeout[0,0])
                de = datetime.datetime.utcfromtimestamp(timeout[-1,-1])
                self.titleString= " %d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (ds.month,ds.day,ds.year,ds.hour+ds.minute/60.0+ds.second/3600.0,\
                    de.month,de.day,de.year,de.hour+de.minute/60.0+de.second/3600.0)
                                     
            try:
                self.pppE
                self.covarE
            except:
                Bmed=5e-5
                self.pppE=(np.array(self.ppp).copy()).tolist(); self.pppE[0]*=Bmed; self.pppE[2]*=Bmed; self.pppE[3]*=Bmed
                self.covarE=(np.array(self.covar).copy()).tolist(); self.covarE[0]*=Bmed*Bmed; self.covarE[1]*=Bmed*Bmed; self.covarE[2]*=Bmed*Bmed
            
            
            if (timeout[-1,-1]-timeout[0,0])/3600.0>36.0:
                if createAll:      
                    self.createPlots(ofname)
                else:
                    self.createPlots_byTime(ofname,txMax=txMax,binByDay=binByDay)                  
            else:
                self.createPlots(ofname)
    
        return
        
    def createPlots_byTime(self,outputFname,txMax=1.0e6,saveout=1,binByDay=0):
         
        vvels1 = self.VectorVels['Vest']; dvvels1 = self.VectorVels['errVest']
        Vmag = self.VectorVels['Vmag']; dVmag = self.VectorVels['errVmag']
        Vdir = self.VectorVels['Vdir']; dVdir = self.VectorVels['errVdir']
 
        evec1 = self.VectorVels['Eest']; devec1 = self.VectorVels['errEest']
        Emag = self.VectorVels['Emag']; dEmag = self.VectorVels['errEmag']
        Edir = self.VectorVels['Edir']; dEdir = self.VectorVels['errEdir']
        
        if self.byGeo==2:
            label = 'Latitude (Deg.)'
            lat = self.VectorVels['Latitude']
            print lat
        else:
            label = 'Magnetic Lat. (Deg.)'
            lat = self.VectorVels['MagneticLatitude']           
               
        MLTtime1 = self.Time['MagneticLocalTime']
        timeout = self.Time['UnixTime'] 
        
        dates=[]; yrs=[]; mths=[]; days=[]; dtime=[]
        for i in range(timeout.shape[0]):
            dates.append([datetime.datetime.utcfromtimestamp(timeout[i,0]),datetime.datetime.utcfromtimestamp(timeout[i,1])])
            yrs.append([dates[i][0].year,dates[i][1].year])
            mths.append([dates[i][0].month,dates[i][1].month])
            days.append([dates[i][0].day,dates[i][1].day])
            dtime.append([dates[i][0].hour+dates[i][0].minute/60.0+dates[i][0].second/3600.0,dates[i][1].hour+dates[i][1].minute/60.0+dates[i][1].second/3600.0])
        yrs=np.array(yrs); mths=np.array(mths); days=np.array(days); dtime=np.array(dtime);
               
        if binByDay:
            r1=datetime.date(dates[0][0].year,dates[0][0].month,dates[0][0].day)
            r2=datetime.date(dates[-1][-1].year,dates[-1][-1].month,dates[-1][-1].day)
            
            # Added 1/7/13 - By Steven Chen
            # Fix for plotting when changing years.
            if int(r1.strftime('%j')) > int(r2.strftime('%j')):
                Ntrecs=int(r2.strftime('%j'))+365 - int(r1.strftime('%j'))+1
            else:
                Ntrecs=int(r2.strftime('%j'))-int(r1.strftime('%j'))+1
            tplus='byDay'
        else:
            Ttotal=(timeout[-1,-1]-timeout[0,0])/3600.0
            Ntrecs=np.ceil(Ttotal/txMax)
            tplus='by' + str(txMax) + 'hr'
                   
        iStart=0; 
        for iTime in range(Ntrecs):
            if iStart>=timeout.shape[0]:
                break;
                        
            if binByDay:
                iEnd=np.where((yrs[:,0]==yrs[iStart,0]) & (mths[:,0]==mths[iStart,0]) & (days[:,0]==days[iStart,0]))[0]
                iEnd=iEnd[-1]
            else:
                iEnd=np.where(timeout[:,-1]<=(timeout[iStart,0]+txMax*3600.0))[0]
                iEnd=iEnd[-1]
            tlim=[iStart,iEnd]
            iStart=iEnd+1
            if Ntrecs>1:
                txtra='-' + tplus + '-' + str(iTime)
            else:
                txtra=''
                
            title= " %d-%d-%d %.3f UT - %d-%d-%d %.3f UT" % (dates[tlim[0]][0].month,dates[tlim[0]][0].day,dates[tlim[0]][0].year,dtime[tlim[0]][0],dates[tlim[-1]][1].month,dates[tlim[-1]][1].day,dates[tlim[-1]][1].year,dtime[tlim[-1],1])
            print title
                
            # velocity vector
            vecPlot = plotVvels.vvelsPlot()
            tvv = vvels1[tlim[0]:(tlim[1]+1)]; dtvv = dvvels1[tlim[0]:(tlim[1]+1)]
            vecPlot.makePlot(timeout[tlim[0]:(tlim[1]+1)],np.mean(MLTtime1[tlim[0]:(tlim[1]+1)],axis=1),lat,tvv[:,:,1],tvv[:,:,0],dtvv[:,:,1],dtvv[:,:,0],\
                title='Vector Velocities' + title,p=self.ppp,sc=self.sc,cax=self.clim,ncols=3,\
                vz=tvv[:,:,2],dvz=dtvv[:,:,2],vzsc=10.0,label=label,geo=self.byGeo)
            if self.saveout:
                print outputFname+'-vvec' +txtra +'.png'
                vecPlot.figg.savefig(outputFname+'-vvec' +txtra +'.png')  
                vecPlot.close()
                
            # velocity magnitude
            magPlot = plotVvels.vvelsMagPlot()            
            magPlot.makePlot(timeout[tlim[0]:(tlim[1]+1)],np.mean(MLTtime1[tlim[0]:(tlim[1]+1)],axis=1),lat,Vmag[tlim[0]:(tlim[1]+1)],Vdir[tlim[0]:(tlim[1]+1)],dVmag[tlim[0]:(tlim[1]+1)],dVdir[tlim[0]:(tlim[1]+1)],
                    title='Vector Velocities' + title,cax1=[0.0,self.clim[1]],cax2=[-180.0,180.0],label=label)
            if self.saveout:
                magPlot.figg.savefig(outputFname+'-vmag' +txtra +'.png')
                magPlot.close()
                
            # electric field vector
            vecPlot = plotVvels.vvelsPlot()
            tvv = evec1[tlim[0]:(tlim[1]+1)]; dtvv = devec1[tlim[0]:(tlim[1]+1)]
            vecPlot.makePlot(timeout[tlim[0]:(tlim[1]+1)],np.mean(MLTtime1[tlim[0]:(tlim[1]+1)],axis=1),lat,tvv[:,:,1]*1e3,tvv[:,:,0]*1e3,dtvv[:,:,1]*1e3,dtvv[:,:,0]*1e3,\
                title='Electric Fields' + title,p=[self.pppE[0]*1e3,self.pppE[1],self.pppE[2]*1e3,self.pppE[3]*1e3],sc=self.sc*5e-2,\
                cax=[self.clim[0]*5e-2,self.clim[1]*5e-2],ncols=3,\
                vz=tvv[:,:,2]*1e3,dvz=dtvv[:,:,2]*1e3,vzsc=10.0,label=label,geo=self.byGeo,units='mV/m',parm='E')
            if self.saveout:
                vecPlot.figg.savefig(outputFname+'-evec' +txtra +'.png')        
                magPlot.close()
            
            # electric field magnitude
            magPlot = plotVvels.vvelsMagPlot()
            magPlot.makePlot(timeout[tlim[0]:(tlim[1]+1)],np.mean(MLTtime1[tlim[0]:(tlim[1]+1)],axis=1),lat,Emag[tlim[0]:(tlim[1]+1)]*1e3,Edir[tlim[0]:(tlim[1]+1)],dEmag[tlim[0]:(tlim[1]+1)]*1e3,dEdir[tlim[0]:(tlim[1]+1)],
                    title='Electric Fields' + title,cax1=[0.0,self.clim[1]*5e-2],cax2=[-180.0,180.0],label=label,units='mV/m',parm='E')
            if self.saveout:
                magPlot.figg.savefig(outputFname+'-emag' +txtra +'.png')  
                magPlot.close()
                
        return
                   
    def createPlots(self,outputFname):
    
        vvels1 = self.VectorVels['Vest']; dvvels1 = self.VectorVels['errVest']
        Vmag = self.VectorVels['Vmag']; dVmag = self.VectorVels['errVmag']
        Vdir = self.VectorVels['Vdir']; dVdir = self.VectorVels['errVdir']
 
        evec1 = self.VectorVels['Eest']; devec1 = self.VectorVels['errEest']
        Emag = self.VectorVels['Emag']; dEmag = self.VectorVels['errEmag']
        Edir = self.VectorVels['Edir']; dEdir = self.VectorVels['errEdir']
        
        if self.byGeo==2:
            label = 'Latitude (Deg.)'
            lat = self.VectorVels['Latitude']
        else:
            label = 'Magnetic Lat. (Deg.)'
            lat = self.VectorVels['MagneticLatitude']
        
        timeout = self.Time['UnixTime']
        MLTtime1 = self.Time['MagneticLocalTime']

        # velocity vector
        vecPlot = plotVvels.vvelsPlot()
        vecPlot.makePlot(timeout,np.mean(MLTtime1,axis=1),lat,vvels1[:,:,1],vvels1[:,:,0],dvvels1[:,:,1],dvvels1[:,:,0],\
            title='Vector Velocities' + self.titleString,p=self.ppp,sc=self.sc,cax=self.clim,ncols=3,\
            vz=vvels1[:,:,2],dvz=dvvels1[:,:,2],vzsc=10.0,label=label,geo=self.byGeo)
        if self.saveout:
            vecPlot.figg.savefig(outputFname+'-vvec.png')  
            #vecPlot.figg.savefig(outputFname+'-vvec.eps')  
            vecPlot.close()      
        
        # velocity magnitude
        magPlot = plotVvels.vvelsMagPlot()
        magPlot.makePlot(timeout,np.mean(MLTtime1,axis=1),lat,Vmag,Vdir,dVmag,dVdir,
                title='Vector Velocities' + self.titleString,cax1=[0.0,self.clim[1]],cax2=[-180.0,180.0],label=label)
        if self.saveout:
            magPlot.figg.savefig(outputFname+'-vmag.png')
            magPlot.close()
            
        # electric field vector
        vecPlot = plotVvels.vvelsPlot()
        vecPlot.makePlot(timeout,np.mean(MLTtime1,axis=1),lat,evec1[:,:,1]*1e3,evec1[:,:,0]*1e3,devec1[:,:,1]*1e3,devec1[:,:,0]*1e3,\
            title='Electric Fields' + self.titleString,p=[self.pppE[0]*1e3,self.pppE[1],self.pppE[2]*1e3,self.pppE[3]*1e3],sc=self.sc*5e-2,\
            cax=[self.clim[0]*5e-2,self.clim[1]*5e-2],ncols=3,\
            vz=evec1[:,:,2]*1e3,dvz=devec1[:,:,2]*1e3,vzsc=10.0,label=label,geo=self.byGeo,units='mV/m',parm='E')
        if self.saveout:
            vecPlot.figg.savefig(outputFname+'-evec.png')        
            vecPlot.close()
        
        # electric field magnitude
        magPlot = plotVvels.vvelsMagPlot()
        magPlot.makePlot(timeout,np.mean(MLTtime1,axis=1),lat,Emag*1e3,Edir,dEmag*1e3,dEdir,
                title='Electric Fields' + self.titleString,cax1=[0.0,self.clim[1]*5e-2],cax2=[-180.0,180.0],label=label,units='mV/m',parm='E')
        if self.saveout:
            magPlot.figg.savefig(outputFname+'-emag.png')            
            magPlot.close()

        return
        
    def setOutPutNames(self,extraApp=''):
        self.outputFnames=[]  
        for Time2Integrate in self.Time2Integrate:
            ofname='%s-%dsec%s' % (self.outputFname,Time2Integrate,extraApp)    
            self.outputFnames.append(ofname)
        return
                          
    def dovels(self):

        # read data
        dat1 = self.h5file.readWholeh5file()

        # output latitude
        x=np.arange(self.plats[0][0],self.plats[0][1],self.plats[0][3])[:,np.newaxis]
        PLAT_OUT=np.concatenate((x,x+self.plats[0][2]),axis=1)
        if len(self.plats[1])>0:
            x=np.arange(self.plats[1][0],self.plats[1][1],self.plats[1][3])[:,np.newaxis]
            x=np.concatenate((x,x+self.plats[1][2]),axis=1)
            PLAT_OUT=np.concatenate((PLAT_OUT,x),axis=0)

        # beamcodes
        BeamCodes=dat1['/']['BeamCodes']

        # only use beams that we were asked to use, if empty,
        # use all beams.
        num_beamcodes = BeamCodes[:,0].size
        if len(self.beams2use) > 0:
            bm_inds = list()
            available_bms = list(BeamCodes[:,0])
            for beamcode in self.beams2use:
                if beamcode in available_bms:
                    bm_inds.append(available_bms.index(beamcode))
            bm_inds = np.array(bm_inds)
        else:
            bm_inds = np.arange(0,num_beamcodes)

        BeamCodes = BeamCodes[bm_inds,:]

        if self.CorrectVap:
            IupB=np.where(BeamCodes[:,0]==upBcode)[0]
            InotUpB=np.where(BeamCodes[:,0]!=upBcode)[0]

        # geomag
        if self.byGeo>0: # geographic
            kpn=dat1['/Geomag']['kn'][bm_inds,:]
            kpe=dat1['/Geomag']['ke'][bm_inds,:]
            kpar=dat1['/Geomag']['kz'][bm_inds,:]
            if self.byGeo==2:
                plat=dat1['/Geomag']['Latitude'][bm_inds,:]
                plong=dat1['/Geomag']['Longitude'][bm_inds,:]
            else:
                plat=dat1['/Geomag']['MagneticLatitude'][bm_inds,:]
                plong=dat1['/Geomag']['MagneticLongitude'][bm_inds,:]
        else: # geomagnetic
            kpn=dat1['/Geomag']['kpn'][bm_inds,:]
            kpe=dat1['/Geomag']['kpe'][bm_inds,:]
            kpar=dat1['/Geomag']['kpar'][bm_inds,:]
            try:
                plat=dat1['/Geomag']['MagneticLatitude'][bm_inds,:]
                plong=dat1['/Geomag']['MagneticLongitude'][bm_inds,:]
            except:
                plat=dat1['/Geomag']['plat']    #???
                plong=dat1['/Geomag']['plong']  #???
        Bn=dat1['/Geomag']['Bx'][bm_inds,:]
        Be=dat1['/Geomag']['By'][bm_inds,:]
        Bz=dat1['/Geomag']['Bz'][bm_inds,:]
        Babs=dat1['/Geomag']['Babs'][bm_inds,:]
        Bmed = np.nanmedian(Babs)
        for i in range(Bmed.ndim):
            Bmed=np.nanmedian(Bmed)		
        self.pppE=(np.array(self.ppp).copy()).tolist(); self.pppE[0]*=Bmed; self.pppE[2]*=Bmed; self.pppE[3]*=Bmed
        self.covarE=(np.array(self.covar).copy()).tolist(); self.covarE[0]*=Bmed*Bmed; self.covarE[1]*=Bmed*Bmed; self.covarE[2]*=Bmed*Bmed

        # fitted params
        ht=dat1['/FittedParams']['Altitude'][bm_inds,:]
        vlos1=dat1['/FittedParams']['Fits'][:,bm_inds,:,0,3] + self.chirp
        vlos1 = np.swapaxes(vlos1,0,1)
        dvlos1=dat1['/FittedParams']['Errors'][:,bm_inds,:,0,3]
        dvlos1 = np.swapaxes(dvlos1,0,1)
        ne1=dat1['/FittedParams']['Ne'][:,bm_inds,:]
        
        # time
        time1=dat1['/Time']['UnixTime']
        doy1=dat1['/Time']['doy']
        dtime1=dat1['/Time']['dtime']+(doy1-doy1[0,0])*24.0
        MLT=dat1['/Time']['MagneticLocalTimeSite']
        yr=dat1['/Time']['Year']
        mon=dat1['/Time']['Month']
        day=dat1['/Time']['Day']
        
        # low densities
        I=np.where((ne1 < self.neMin))
        vlos1[I]=np.nan
        dvlos1[I]=np.nan
    
        # just do a portion of data
        if len(self.zoomWhole)!=0:
            I=np.where((dtime1[:,0]>=self.zoomWhole[0]) & (dtime1[:,1]<=self.zoomWhole[1]))[0]	
            vlos1=vlos1[I]
            dvlos1=dvlos1[I]
            time1=time1[I]
            doy1=doy1[I]		
            dtime1=dtime1[I]		
            MLT=MLT[I]
            yr=yr[I]
            mon=mon[I]
            day=day[I]
    
        # title str
        yr=yr[0,0]; mon=mon[0,0]; day=day[0,0]
        self.titleString=' %d-%d-%d' % (mon, day, yr)
    
        (Nbeams,Nhts)=ht.shape

        if self.CorrectVap:
            Nbeams=Nbeams-1
            # up B
            vlos1upB=vlos1[:,IupB[0],:]
            dvlos1upB=dvlos1[:,IupB[0],:]
            htupB=ht[IupB[0],:]
            # not
            vlos1=vlos1[:,InotUpB,:]
            dvlos1=dvlos1[:,InotUpB,:]
            ht=ht[InotUpB,:]	
            kpn=kpn[InotUpB,:]
            kpe=kpe[InotUpB,:]
            kpar=kpar[InotUpB,:]
            plat=plat[InotUpB,:]
            plong=plong[InotUpB,:]    
                  
        for itimeInt in range(len(self.Time2Integrate)):
            
            Time2Integrate = self.Time2Integrate[itimeInt]
            ofname=self.outputFnames[itimeInt]       
            
            # output arrays
            timeout=[]
            dtimeout=[]
            MLTtime1=[]
            
            done=0 # flag to say when we are done
            Irec=0 # record counter
            IIrec=0 # record counter

            while not done:
                Irecs=np.where((time1[:,0]>=time1[Irec,0]) & (time1[:,0]<=(time1[Irec,0]+Time2Integrate)))[0]

                tvlos=vlos1[Irecs,:,:]
                tdvlos=dvlos1[Irecs,:,:]
                        
                if self.CorrectVap:
                    vUpB=vlos1upB[Irecs]
                    dvUpB=dvlos1upB[Irecs]
                    for i in range(Nbeams):
                        tupB=np.interpolate.interp1d(np.squeeze(htupB),vUpB,bounds_error=0,fill_value=0.0)(ht[i,:])
                        tvlos[:,i,:]=tvlos[:,i,:]-tupB*kpar[i,:]
                        tdvlos[:,i,:]=np.sqrt(np.power(tdvlos[:,i,:],2.0)+np.power(dvUpB*kpar[i,:],2.0))
                
                # line of sight velocities and errors
                vlosin=np.reshape(np.squeeze(tvlos),(Nhts*Nbeams*len(Irecs)))
                dvlosin=np.reshape(np.squeeze(tdvlos),(Nhts*Nbeams*len(Irecs)))

                # k vector (this is in mag coords)
                kin=np.zeros((Nbeams,Nhts,3),dtype=kpn.dtype)
                kin[:,:,0]=kpn
                kin[:,:,1]=kpe
                kin[:,:,2]=kpar
                kin=np.reshape(np.repeat(kin[np.newaxis,:,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams,3))

                ekin=np.zeros((Nbeams,Nhts,3),dtype=kpn.dtype)
                ekin[:,:,0]=(-Be*kpar-Bz*kpe)/Babs**2.0
                ekin[:,:,1]=(Bz*kpn+Bn*kpar)/Babs**2.0
                ekin[:,:,2]=(Bn*kpe-Be*kpn)/Babs**2.0
                ekin=np.reshape(np.repeat(ekin[np.newaxis,:,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams,3))
                
                # magnetic field
                bin=np.reshape(np.repeat(Babs[np.newaxis,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams))
                
                # mag lat and long, altitude
                platin=np.reshape(np.repeat(plat[np.newaxis,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams))
                plongin=np.reshape(np.repeat(plong[np.newaxis,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams))
                htin=np.reshape(np.repeat(ht[np.newaxis,:,:],len(Irecs),axis=0),(len(Irecs)*Nhts*Nbeams))

                # compute vectors
                (plat_out1,Vest,dVest,vVestAll,Nmeas,vchi2)=vvels.compute_velvec2(PLAT_OUT,vlosin,dvlosin,kin,platin,plongin,htin,htmin=self.minAlt*1000,htmax=self.maxAlt*1000,covar=self.covar,p=self.ppp)
                (plat_out1,Eest,dEest,vEestAll,Nmeas1,echi2)=vvels.compute_velvec2(PLAT_OUT,vlosin,dvlosin,ekin,platin,plongin,htin,htmin=self.minAlt*1000,htmax=self.maxAlt*1000,covar=self.covarE,p=self.ppp)
                Eest[:,-1]*=-1

                timeout.append([time1[Irecs[0],0],time1[Irecs[-1],1]])
                dtimeout.append([dtime1[Irecs[0],0],dtime1[Irecs[-1],1]])
                MLTtime1.append([MLT[Irecs[0],0],MLT[Irecs[-1],1]])
                if IIrec>0:
                    while MLTtime1[IIrec][0]<MLTtime1[IIrec-1][0]:
                        MLTtime1[IIrec][0]=MLTtime1[IIrec][0]+24.0
                        MLTtime1[IIrec][1]=MLTtime1[IIrec][1]+24.0

                if IIrec==0:
                    vvels1=Vest[np.newaxis,:,:]
                    dvvels1=dVest[np.newaxis,:,:]
                    varvvels1=vVestAll[np.newaxis,:,:,:]
                    Nall1=Nmeas[np.newaxis,:]
                    chi2 = vchi2[np.newaxis,:]
                    evec1=Eest[np.newaxis,:,:]
                    devec1=dEest[np.newaxis,:,:]
                    varevec1=vEestAll[np.newaxis,:,:,:]
                else:
                    vvels1=np.concatenate((vvels1,Vest[np.newaxis,:,:]),axis=0)
                    dvvels1=np.concatenate((dvvels1,dVest[np.newaxis,:,:]),axis=0)
                    varvvels1=np.concatenate((varvvels1,vVestAll[np.newaxis,:,:,:]),axis=0)
                    Nall1=np.concatenate((Nall1,Nmeas[np.newaxis,:]),axis=0)
                    chi2 = np.concatenate((chi2,vchi2[np.newaxis,:]),axis=0)
                    evec1=np.concatenate((evec1,Eest[np.newaxis,:,:]),axis=0)
                    devec1=np.concatenate((devec1,dEest[np.newaxis,:,:]),axis=0)
                    varevec1=np.concatenate((varevec1,vEestAll[np.newaxis,:,:,:]),axis=0)

                Irec=Irecs[-1]+1 # increment counters
                IIrec=IIrec+1

                if Irec>=time1.shape[0]:
                    done=1



            MLTtime1=np.array(MLTtime1)
            timeout=np.array(timeout)
            dtimeout=np.array(dtimeout)
            
            # magnitude and direction
            Vmag = np.sqrt( np.power(vvels1[:,:,0],2.0) + np.power(vvels1[:,:,1],2.0) ).real
            dVmag = np.sqrt( np.power(dvvels1[:,:,0],2.0)*np.power(vvels1[:,:,0]/Vmag,2.0) + np.power(dvvels1[:,:,1],2.0)*np.power(vvels1[:,:,1]/Vmag,2.0) ).real
            Vdir = 180.0/np.pi*np.arctan2(vvels1[:,:,1],vvels1[:,:,0]).real
            dVdir = 180.0/np.pi*((1.0/np.absolute(vvels1[:,:,0]))*(1.0/(1.0+np.power(vvels1[:,:,1]/vvels1[:,:,0],2.0)))*np.sqrt(np.power(dvvels1[:,:,1],2.0)+np.power(vvels1[:,:,1]/vvels1[:,:,0]*dvvels1[:,:,0],2.0))).real		
            
            Emag = np.sqrt( np.power(evec1[:,:,0],2.0) + np.power(evec1[:,:,1],2.0) ).real
            dEmag = np.sqrt( np.power(devec1[:,:,0],2.0)*np.power(evec1[:,:,0]/Emag,2.0) + np.power(devec1[:,:,1],2.0)*np.power(evec1[:,:,1]/Emag,2.0) ).real
            Edir = 180.0/np.pi*np.arctan2(evec1[:,:,1],evec1[:,:,0]).real
            dEdir = 180.0/np.pi*((1.0/np.absolute(evec1[:,:,0]))*(1.0/(1.0+np.power(evec1[:,:,1]/evec1[:,:,0],2.0)))*np.sqrt(np.power(devec1[:,:,1],2.0)+np.power(evec1[:,:,1]/evec1[:,:,0]*devec1[:,:,0],2.0))).real		
            
            ### Set up output dicts
            self.setOutDicts(dat1)
            self.Params['IntegrationTime']=Time2Integrate
            # Time
            self.Time['UnixTime'] = timeout
            self.Time['dtime'] = dtimeout
            self.Time['MagneticLocalTime'] = MLTtime1
            # Vector Vels
            if self.byGeo==2:
                self.VectorVels['Latitude']=PLAT_OUT
                self.VectorVels['MagneticLatitude']=[]
            else:
                self.VectorVels['Latitude']=[]
                self.VectorVels['MagneticLatitude']=PLAT_OUT
            self.VectorVels['Nmeas'] = Nall1.astype('int32')
            self.VectorVels['Chi2'] = chi2
            self.VectorVels['Vest'] = vvels1; self.VectorVels['errVest'] = dvvels1
            self.VectorVels['varVest'] = varvvels1
            self.VectorVels['Vmag'] = Vmag; self.VectorVels['errVmag'] = dVmag
            self.VectorVels['Vdir'] = Vdir; self.VectorVels['errVdir'] = dVdir
            self.VectorVels['Eest'] = evec1; self.VectorVels['errEest'] = devec1
            self.VectorVels['varEest'] = varevec1
            self.VectorVels['Emag'] = Emag; self.VectorVels['errEmag'] = dEmag
            self.VectorVels['Edir'] = Edir; self.VectorVels['errEdir'] = dEdir
            
            ### Output file
            if self.saveout:
                self.createOutputFile(ofname)
                
            ### Plot output
            if self.makeplot:            
                if (timeout[-1,-1]-timeout[0,0])/3600.0>36.0:
                    self.createPlots_byTime(ofname,binByDay=1)
                else:
                    self.createPlots(ofname)                
        
        return
        
