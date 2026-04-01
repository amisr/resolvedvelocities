import numpy
import datetime
import os
import pickle
import tables
import glob
import urllib
import requests
import shutil
import time

class GeophysicalParameters:

    def __init__(self, DirIn, InFile, MagFileLocation, DelayTime = 2):
        # self.DirIn = '/Users/srkaeppler/research/data/NSF_ERegionNeutralWind/geophys_params/'
        # self.InFile = 'KP_AP_F107_2006_2017.pkl'

        self.DirIn = DirIn
        self.InFile = InFile
        self.MagFileLocation = MagFileLocation
        self.DelayTime = DelayTime

        with open(os.path.join(self.DirIn,self.InFile), 'rb') as f:
            self.AP_KP_data = pickle.load(f, encoding="latin1")

        self.AP_KP_UnixStartTime = numpy.array([x['UnixStartTime'] for x in self.AP_KP_data])
        self.AP_KP_UnixEndTime = numpy.array([x['UnixEndTime'] for x in self.AP_KP_data])
        return

    def GetGeophysicalParameters(self,tUnixIn):
        """
        This program will get the AP, KP, AE, f107 and produce an output List
        Variables can then be assigned to the proper output format

        Inputs: unixtime arry
        Output: numpy dictionary of AP, KP, AE, F107, etc

        """

        N = tUnixIn.shape[0]
        KP = numpy.zeros(N)*numpy.nan
        AP = numpy.zeros(N)*numpy.nan
        AE = numpy.zeros(N)*numpy.nan
        AL = numpy.zeros(N)*numpy.nan
        AU = numpy.zeros(N)*numpy.nan
        SymH = numpy.zeros(N)*numpy.nan
        APmean = numpy.zeros(N)*numpy.nan
        KPsum = numpy.zeros(N)*numpy.nan
        F107 = numpy.zeros(N)*numpy.nan

        # get the correct year for the OMNI data
        yearStart = datetime.datetime.utcfromtimestamp(tUnixIn[0]).year
        yearEnd = datetime.datetime.utcfromtimestamp(tUnixIn[-1]).year

        if yearStart == yearEnd:
            fileAE  = glob.glob(os.path.join(self.DirIn, 'OMNI_5MIN_%s_04232020.h5'%yearStart))
            # print fileAE
            if len(fileAE) == 1:
                self.AE_SymH_data = tables.open_file(fileAE[0])
                # found a comment that argmin is very fast with numpy arrays - most time
                # spent converting data:
                #https://stackoverflow.com/questions/26816590/faster-alternatives-to-numpy-argmax-argmin-which-is-slow
                AE_tunix = numpy.array(self.AE_SymH_data.root.UnixTime.read())


        for itime in range(len(tUnixIn)):
            qAP = numpy.where((self.AP_KP_UnixStartTime < tUnixIn[itime]) & \
                              (self.AP_KP_UnixEndTime >= tUnixIn[itime]) )[0][0]

            if qAP:
                KP[itime] = self.AP_KP_data[qAP]['KP']
                AP[itime] = self.AP_KP_data[qAP]['AP']
                KPsum[itime] = self.AP_KP_data[qAP]['KPsum']
                APmean[itime] = self.AP_KP_data[qAP]['APmean']
                F107[itime] = self.AP_KP_data[qAP]['f107']

            # now grab the AE parameters
            # this is a slow                     )
            qAE = numpy.argmin(numpy.abs(AE_tunix-tUnixIn[itime]))
            if qAE:
                AE[itime] = self.AE_SymH_data.root.AE[qAE]
                AL[itime] = self.AE_SymH_data.root.AL[qAE]
                AU[itime] = self.AE_SymH_data.root.AU[qAE]
                SymH[itime] = self.AE_SymH_data.root.SymH[qAE]

        outDict = {}
        outDict['KP'] = KP
        outDict['AP'] = AP
        outDict['KPsum'] = KPsum
        outDict['APmean'] = APmean
        outDict['F107'] = F107
        outDict['AE'] = AE
        outDict['AL'] = AL
        outDict['AU'] = AU
        outDict['SymH'] = SymH

        self.AE_SymH_data.close()
        return outDict


    def ProcessGroundMag(self,tUnixStart,tUnixEnd):

        """
        Processing code to process local magnetometer data as another proxy
        Original code is located in the ../geophy_Params folder

        I decided to take the beginning time and the end time of an experiments
        and grab all of the data in between there.
        From that, I also decimate the data by about a factor of 10, so it has
        about a 10 second time cadance which should be sufficient
        """

        tStart = datetime.datetime.utcfromtimestamp(tUnixStart)
        tEnd = datetime.datetime.utcfromtimestamp(tUnixEnd)
        t1970 = datetime.datetime(1970,1,1,00,00,00)
        # rootStr = 'https://www.asf.alaska.edu/magnetometer/archive/'
        # updated on 4/10/2020
        rootStr = 'https://www.gi.alaska.edu/api/magnetometer/DATA/www/'

        datetimeList = numpy.unique(numpy.array([tStart,tEnd]))
        UnixTime = []
        H = []
        D = []
        Z = []
        urlList = []
        for idt in datetimeList:
            tmpUnixTime = (datetime.datetime(idt.year,idt.month,idt.day,00,00,00)-t1970).total_seconds()
            iurl = os.path.join(rootStr,'%04d'%idt.year,\
                                '%02d'%idt.month,'%02d'%idt.day,\
                                'poker','poker_%04d_%02d_%02d.csv'%(idt.year,idt.month,idt.day))
            urlList.append(iurl)

        urlList = numpy.unique(numpy.array(urlList))
        for iurl in urlList:
            print(iurl)
            # check if file exists on remote server, if not skip

            # had an idea, will check if file exists locally first.
            # if it does, use that, if not, request it and save a copy to a directory

            localFileLocation = os.path.join(self.MagFileLocation,iurl.split('/')[-1])

            if os.path.isfile(localFileLocation):
                print('Grabbing Mag File Locally')
                dataIn = numpy.loadtxt(localFileLocation, delimiter=',', usecols=[1,2,3,4])
                UnixTime.extend(tmpUnixTime + numpy.array(dataIn[::30,0])*3600.)
                H.extend(dataIn[::30,1])
                D.extend(dataIn[::30,2])
                Z.extend(dataIn[::30,3])
            else:
                print('Grabbing Mag File From Server')
                time.sleep(self.DelayTime)
                mkrequest = requests.head(iurl)
                if mkrequest.status_code == requests.codes.ok:
                    inFile =urllib.urlretrieve(iurl)
                    dataIn = numpy.loadtxt(inFile[0], delimiter=',', usecols=[1,2,3,4])
                    # copy the file to the local location
                    shutil.copy2(inFile[0],localFileLocation)
                    os.remove(inFile[0])
                    del inFile

                    UnixTime.extend(tmpUnixTime + numpy.array(dataIn[::10,0])*3600.)
                    H.extend(dataIn[::10,1])
                    D.extend(dataIn[::10,2])
                    Z.extend(dataIn[::10,3])


        UnixTime = numpy.array(UnixTime)
        H = numpy.array(H)
        D = numpy.array(D)
        Z = numpy.array(Z)

        qsort = numpy.argsort(numpy.unique(UnixTime))
        UnixTime = UnixTime[qsort]
        H = H[qsort]
        D = D[qsort]
        Z = Z[qsort]

        outDict = {}
        outDict['GroundMag_H'] = H
        outDict['GroundMag_D'] = D
        outDict['GroundMag_Z'] = Z
        outDict['GroundMag_UnixTime'] = UnixTime

        return outDict

if __name__ == '__main__':
    # test script

    t1970 = datetime.datetime(1970,1,1,0,0,0)
    t2014 = datetime.datetime(2014,1,1,00,00,1)

    t = (t2014-t1970).total_seconds()+numpy.arange(0,3600*36,180)
    GeoPhyClass = GeophysicalParameters()
    outDict = GeoPhyClass.GetGeophysicalParameters(t)
    magDict = GeoPhyClass.ProcessGroundMag(t[0],t[-1])
