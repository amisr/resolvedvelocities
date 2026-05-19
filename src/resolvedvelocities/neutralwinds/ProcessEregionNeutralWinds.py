"""
Process Eregion Neutral Winds code.  Developed in support of PFISR NSF E-region
Neutral winds project in collaboration with Clemson University.

Author: Stephen Kaeppler
Email: steve.kaeppler@gmail.com

Purpose: This is the overarching program and functions which process the
E region neutral winds from the fitted AC and LP data from PFISR.
This is a conversion fo process_eregwinds_srk.py which was originally written by
Nicolls into a more formal python class structure.

2017-10-05 - v0.2

The ProcessEregionNeutralWinds.py file has been validated against process_eregwinds_srk.py
using 20161121.001_ac_3min-fitcal.h5, 20170301.013_ac_3min-fitcal.h5, 20170302.001_ac_3min-fitcal.h5.
The program to run these is ComparePrograms.py.  At this point these  program match.
I am going to start diverging the code base, first subtly in the Joule Heating
since I found that Mike just looped over Nbeams, which isn't quite right, you need to loop
over the beams that were selected.

Changes from this point forward will produce different results.

2017-10-10 - v0.3.2017.10.10

Version v0.3, I made some IO changes but I may start processing some data with this version.

Version v0.4 - lots of small edits made to the IO and the plotting software.  It all seems to work
I have also included the SNR and Ne into the monthly plots and other information.
Made processing smoother.

03 13 2018 - added solar local time converion

v0.4.1 - 09 08 2018 added some ability to extract out the raw electron and SNR densities for each altitude bin
v0.4.2 - 10 15 2018 added in obtaining the F-region flows - want to check against the electric field.
v0.4.3 - 10 29 2018 added in some more altitude into the Joule Heating so I can make better figures
v0.4.4 - 11 20 2018 made some pretty major changes to IO to include consistent calculation of
Pedersen conductivity from FastConductivity.py.  Made some changes to the Joule heating calculation and checked
formulas.  It is worth checking again.

v0.4.5 - 11 20 2018: added in Hall and Pedersen conductivities from fitted electron density data.
v0.4.6 - 12 03 2018: Tried to fix some of the double counting and time problems in testMakeMonthlyh5

05 22 2019: added some statements to bypass the geophysical parameters.  Also need in config file now.
            Additionally wrote in IOEregionwinds a try except statement

07 29 2019: Running the code for the 06 data reprocessed by Ashton

v0.5.0 - 10-15-2019: put in some filtering on the LOS velocity discharging bad Chi square and bad error codes on the fit.

v0.5.1 - 10-23-2019: changed chi square to 0.01 for lower boundary

v0.5.3 - 12-02-2019: Added in that now passing in the Chi2 and Fitcode filtering by Config file
            Bigger change that I am scaling the AC dVlos by some sort of factor while Ashton figures this out.
            We decided that a conversative scaling would be to reduce the dVLOS by 1/sqrt(10).
            The chi square produced in the data Ashton sent me typically was around 0.01, so the uncertaintiies on the LOS velocities
            may be over estimated.  So we are just changing this as a temporary fix while Ashton fixes the uncertainty estimation.

v0.5.5 02 01 2020 - Added in calculation of Coriolis, Centrifugal, and Lorentz forcing
v0.5.5 02 10 2020 - Added a correction to qvert so that way I can calculate the lorentz term.
                    Found an error where qvert = 0 in the if statement goes to false.

v0.5.5 02 15 2020 - Put in  nuInscaler into the main program, scaling ALL kappas by the scaler number

v0.5.6 02 28 2020 -- Added some more vlos diagostics and the calculation of the scale height. Added Altitude offset

v0.5.6.2020.03.12_nuin_fracoff - testing putting in the Brekke formula for ion neutral collision frequency and took out frac

v0.5.7.2020.04.10 - Put in Ashton's revised ion neutral collision frequency formulas into IO.
                    Also wrote a testscript and at least for the file I used was only different by 2.5%.
                    Revised where the mag data is being pulled from since the URL is deprecated
                    Added in Kappa which is now being interpolate - plan to see where kappa =1 is located for the paper.
                    commented out nuin scaler just so I am not chasing my tail

                    test v0.5.7.2020.04.13_org commented back in original ion neutral collision frequency method
                    possible mistake that not summing up properly.

                    test v0.5.7.2020.04.13_newnuin_orgsum_noTr800 - new formula for nuin except took off Tr>800.
                    I expect this should be almost the same as before since the formulas are basically the same.
                    did the original sum using frac[0] and frac[1] want to see if I am underestimating

                    v0.5.7.2020.04.13_newnuin_orgsum_yesTr800 - same as above except now including Tr>800.

                    'v0.5.7.2020.04.13_newnuin_newsum_noTr800' - using the new sum now and new col freq

                    v0.5.7.2020.04.13_updatedorg - updated original uses original method but including the NO term

v0.6.0.2020.04.15 -- Now think I have the new ion neutral collision frequency working and validated.
                    Found a mistake in how I was calculating the ion neutral collision frequency that
                    the fraction weight I was using only included the O+ and O2+ terms and not NO+
                    Turns out I was basically weighting by about 0.5, so I was effectively reducing the
                    ion neutral collision frequency by about a factor of 0.5 or less...
                    From this point forward need to start using any results from > v0.6
                    This revision has changed previous results signi

v0.6.0.2020.04.21 -- updated to now include the temperature correction for the O2+

v0.6.1.2020.04.23 -- made a number of changes to the geomagnetic files and reprocessed from CDAweb.
                     Wrote new code to be able to process the files from CDAweb in the new format.
                     Also changed the geomagnetic data files

v0.6.1.2020.06.07 -- changed the generation of Monthly files to hopefully be in order now
                    Added in missingIPY files given to me by ashton, maybe improve data covarege
                    Some work going to need to be done to make sure that all of the 10, 15, and 20 minute data are there.

v0.6.1.2020.06.15_Weijia -- Updated the data for Weijia's study in particular since we are missing a lot of IPY data for 02-04 2013 and 2014.

'v0.6.2.2020.07.01' -- Updated the data with new IPY27 mode for 2013 and 2014 Ashton processed.  Also now put in mechanical Joule heating term.
                        Put in the conductance and conductivity now too.

v0.6.2.2020.07.30 -- Made some changes to IO since Weijia noticed the mechanical heating terms were missing from the monthly files.

v0.6.3.2020.10.19 -- Tried to elimated all extra instance of nuinscaler, and also output that variable.  Added in variables
                    To get the Ti, Tn, ion neutral collision frequency along the vertical beam for diagnostic purposes
                    included dVest for F-region plasma drifts for Rafael

v0.6.4.2020.11.20 -- Extracted some more parameters including F107 and the Hall and Pedersen Drags

v0.6.4.2021.07.21 -- Final Run of data for NSF project

v0.6.5.2021.08.06 -- Added in a bunch of the thayer calculations after the reviewer comments we got on Weijia's paper
                    This only affects the Joule heating calculations.

v0.6.6.2022.05.03 -- Added full altitude for alternating code data for Joule heating and Pedersen Conductivity

"""
#import sys
#sys.path.append("./GeophysicalParameters")

from .IOEregionWinds import IOEregionWinds
from .CalculateEregionWinds import CalculateEregionWinds as CalculateEregionWinds
from .GeophysicalParameters import GeophysicalParameters as GeophysicalParameters
import numpy
from scipy.interpolate import interp1d
from scipy.constants import elementary_charge
from .tools.loggerinit.LoggerInit import *
from .tools.configreader.ConfigReader import *
from . import summary_plots
import os
import sys
import datetime
import time
import copy
from argparse import ArgumentParser, RawDescriptionHelpFormatter

try:
    from numpy import trapz
except ImportError:
    from numpy import trapezoid as trapz



VersionNumber = 'v0.6.6.2022.05.03'

config_file_help = """Calculate 3D resolved neutral winds from the LoS measurments
in a fitted AMISR file in altitude bins."""




class ProcessEregionNeutralWinds:

    def __init__(self, configFile):
        # check if config file exists
        # parse the config file
        # set the neutral wind altitude grid

        # LL - replace with standard config reader function?
        self.ConfigReader = ConfigReader()

        # LL - does any of this actually need to be defined in advance??
        self.ProcessingListMaster = ['NW_MinAltitude', 'NW_MaxAltitude','NW_DeltaAltitude',\
                                    'LP_MinAltitude', 'LP_MaxAltitude', 'Elevation_Min', \
                                    'Elevation_Max', 'Azimuth_Min','Azimuth_Max', \
                                    'TestAlgorithm','LocalTimeZone','dVLOSScaler','nuinScaler', 'AltOffset']
        self.VVelsListMaster = ['Velocity_ModelCovariance', 'Velocity_MaxAltitude', \
                                'Velocity_MinAltitude', 'FracErrorOffset', 'FracErrorThreshold', \
                                'AbsoluteErrorThreshold']

        self.InvertWindsListMaster = ['AbsoluteErrorFilter', 'Wind_ModelCovariance', \
                                        'Efield_ModelCovariance']

        self.GeoPhysListMaster = ['ByPass','Path2Data','Path2MagFiles','KP_AP_File','DelayTime']

        #self.Default = ['RootFileDirectory']
        self.Default = ['AC_FileName', 'LP_FileName', 'OutputFileName', 'OutputPath']


        self.DataExclude = ['Chi2Max','Chi2Min','FitCodeMax','FitCodeMin']

        self.configFile = configFile
        if configFile:
            #if os.path.isfile(configFile):
            self.config = self.ConfigReader.read(configFile)
            print(self.config.keys())
            # LL - skip config checking - errors raised regardless and it removes a bunch of code infrastructure to keep up to date?
            self._check_config()
            #else:
            #    raise Exception('Error: Check config file parameters')

            # parse config file at this point

        else:
            # The situation of no config file specified should have been caught several times before?
            raise Exception('Error: No Config File Specified - exiting')
            sys.exit()
            #if no config file specified use standard inputs

        self.nuinScaler = 1.0
        #instantiation of classes
        if 'nuinScaler' in self.config['PROCESSING'].keys():
            self.nuinScaler = self.config['PROCESSING']['nuinScaler']

        self.AltOffset = 0.
        # it is a bad idea to keep this in because everything will get averaged out
        # if 'AltOffset' in self.config['PROCESSING'].keys():
        #     self.AltOffset = self.config['PROCESSING']['AltOffset']

        self.io = IOEregionWinds(nuinScaler=self.nuinScaler)
        self.winds = CalculateEregionWinds()
#        self.geoparams = GeophysicalParameters.GeophysicalParameters(self.config['GEOPHYSICALPARAMETERS']['Path2Data'], \
#                                                                    self.config['GEOPHYSICALPARAMETERS']['KP_AP_File'], \
#                                                                    self.config['GEOPHYSICALPARAMETERS']['Path2MagFiles'], \
#                                                                    self.config['GEOPHYSICALPARAMETERS']['DelayTime'])

        return



    def _check_config(self):
        """
        This is a private function designed to check the config file
        IF something is missing, need to abort the program
        """
        k=0
        for outterkey in self.config.keys():
            for ikey in self.config[outterkey].keys():
                if ikey in self.ProcessingListMaster:
                    # self.ProcessingTruthTable[k] = True
                    print('ikey pass:', ikey)
                elif ikey in self.InvertWindsListMaster:
                    print('ikey pass:', ikey)
                elif ikey in self.VVelsListMaster:
                    print('ikey pass', ikey)
                elif ikey in self.GeoPhysListMaster:
                    print('ikey pass:', ikey)
                elif ikey in self.Default:
                    print('ikey pass:', ikey)
                elif ikey in self.DataExclude:
                    print('ikey pass:', ikey)
                else:
                    raise ValueError('%s -- Required Key not in Config File'%ikey)
                k+=1
        # print self.ProcessingTruthTable

        return





    #def run(self,fname_ac, fname_lp):
    def run(self):

        fname_ac = self.config['DEFAULT']['AC_FileName']
        fname_lp = self.config['DEFAULT']['LP_FileName']

        # make sure that the files exist
        if os.path.isfile(fname_ac):
            acDict= self.io.MakeInputDictionary(fname_ac)
        else:
            raise Exception('Error: Alternating Code Input File does not exist')

        if os.path.isfile(fname_lp):
            lpDict= self.io.MakeInputDictionary(fname_lp)
        else:
            raise Exception('Error: Long Pulse Input File does not exist')


        #oname = fname_ac.split('/')[-2]
        #FnameOut = fname_ac.split('/')[-1][0:-3]
        ## check that the output file exists/make a directory
        #FullPathOut = os.path.join(self.config['DEFAULT']['RootFileDirectory'],VersionNumber,oname)
        FullPathOut = self.config['DEFAULT']['OutputPath']
        FnameOut = self.config['DEFAULT']['OutputFileName']

        if not os.path.exists(FullPathOut):
            os.makedirs(FullPathOut)

        # LL - not needed?
        # make a 'lock file'
        lock_location = os.path.join(FullPathOut,FnameOut+'.lock')
        with open(lock_location, 'w') as f:
            f.write('')

        # SRK - setting an altitude grid
        # Dealing with the altitude grids
        LPminAlt = self.config['PROCESSING']['LP_MinAltitude']
        LPmaxAlt = self.config['PROCESSING']['LP_MaxAltitude']
        tm1=numpy.arange(self.config['PROCESSING']['NW_MinAltitude'], \
                         self.config['PROCESSING']['NW_MaxAltitude'],\
                         self.config['PROCESSING']['NW_DeltaAltitude'])
        #	tm2=scipy.arange(200.0,400.0,50)
        #	htout=scipy.zeros((tm1.size+tm2.size,2),dtype='float64')
        #	htout[:,0]=scipy.concatenate((tm1,tm2))
        #	htout[:,1]=scipy.concatenate((tm1+tm1[1]-tm1[0],tm2+tm2[1]-tm2[0]))
        htout=numpy.zeros((tm1.size,2),dtype='float64')
        htout[:,0]=tm1
        htout[:,1]=tm1+tm1[1]-tm1[0]
        if self.AltOffset:
            htout[:,1] = htout[:,1]+self.AltOffset
        htoutm=numpy.mean(htout,axis=1)

        print('htout', htout)

        # LL - This creates all the empty arrays - would save a lot of space to skip...
        # make the output dictionary
        outDict = self.io.MakeOutputDictionary(acDict,htout)
        #outDict = dict()


        # handling the beams
        BMCODES = acDict['BMCODES']
        Ibm = numpy.where((BMCODES[:,2] >= self.config['PROCESSING']['Elevation_Min']) & \
                            (BMCODES[:,2] <= self.config['PROCESSING']['Elevation_Max']) & \
                            (BMCODES[:,1] >= self.config['PROCESSING']['Azimuth_Min']) & \
                            (BMCODES[:,1] <= self.config['PROCESSING']['Azimuth_Max']))[0]
        Nbeams1 = len(Ibm)
        # print Nbeams1

        # find the vertical beam, if it exists
        qvert = numpy.where(numpy.abs(BMCODES[:,2]-90.) \
                == numpy.min(numpy.abs((BMCODES[:,2]-90.))) )[0]


        # getting time in order
        outDict['UnixTime'] = acDict['UnixTime'].copy()
        outDict['MeanUnixTime'] = numpy.mean(acDict['UnixTime'], axis=1)
        UTDecHrs = numpy.mean(acDict['dtime'], axis=1)
        outDict['UTDecHrs'] = UTDecHrs.copy()
        # do local time, will need to make this generalizable through config
        #03 13 2018 added solar local time conversion
        #original # LocalDecHrs = UTDecHrs + self.config['PROCESSING']['LocalTimeZone']
        #solar local time
        LocalDecHrs = UTDecHrs + acDict['SiteLongitude']/15.
        qtime = numpy.where(LocalDecHrs < 0.)[0]
        LocalDecHrs[qtime] = LocalDecHrs[qtime]+24.
        outDict['LocalDecHrs'] = LocalDecHrs
        outDict['MLTDecHrs'] = numpy.mean(acDict['MLT'], axis=1)

        # a few other odds and ends
        outDict['Altitude'] = htoutm
        #outDict['AltitudeJH'] = htoutm
        statusArr = numpy.ones(acDict['UnixTime'].shape[0])

        ## added on 11 20 2018
        #outDict['PedersenConductivity'] = acDict['PedersenConductivity']
        #outDict['HallConductivity'] = acDict['HallConductivity']
        #outDict['F107A'] = acDict['F107A']
        #outDict['F107Raw'] = acDict['F107Raw']


        ## empty arrays that will be filled
        #Timeby3 = numpy.zeros((acDict['UnixTime'].shape[0],3), dtype='float64')*numpy.nan
        ##Timeby3List = ['Efield', 'errEfield', 'VestGmag_300km','errVestGmag_300km']
        #Timeby3List = ['Efield', 'errEfield']
        #for ii in Timeby3List:
        #    outDict[ii] = copy.copy(Timeby3)

        ## time x altitude x 3 - NW grid
        #TimebyAltitudeby3 = numpy.zeros((acDict['UnixTime'].shape[0],htoutm.shape[0],3),\
        #                    dtype='float64')*numpy.nan
        ##TimebyAltitudeby3List = ['WindGmag', 'errWindGmag', 'WindGeo', 'errWindGeo', \
        ##                        'VestGmag', 'errVestGmag','VestGeo','errVestGeo', \
        ##                        'Coriolis', 'Centrifugal', 'Lorentz']
        #TimebyAltitudeby3List = ['WindGmag', 'errWindGmag', 'WindGeo', 'errWindGeo', \
        #                        'VestGmag', 'errVestGmag','VestGeo','errVestGeo']
        #for ii in TimebyAltitudeby3List:
        #    outDict[ii] = copy.copy(TimebyAltitudeby3)

        #NeShape2 = numpy.zeros((acDict['Ne'].shape[0],acDict['Ne'].shape[1],acDict['Ne'].shape[2],2), \
        #            dtype='float64')*numpy.nan
        #NeShape2List = ['CurrentGmag', 'errCurrentGmag']

        #for ii in NeShape2List:
        #    outDict[ii] = copy.copy(NeShape2)


        for itime in range(acDict['UnixTime'].shape[0]):
            # getting all the Alternating code

            print('itime', itime,acDict['UnixTime'].shape[0])
            print('datetime', datetime.datetime.utcfromtimestamp(acDict['UnixTime'][itime,0]))

            AllBabs = numpy.ravel(acDict['Babs'][Ibm,:])
            Alldip = numpy.ravel(acDict['dip'][Ibm,:])
            Alldec = numpy.ravel(acDict['dec'][Ibm,:])
            Allk = numpy.reshape(acDict['k1'][Ibm,:,:], (Nbeams1*acDict['NAltitudes'],3))
            AllAltitude = numpy.ravel(acDict['Altitude'][Ibm,:])/1000. # convert to km
            AllVlos = numpy.ravel(acDict['Vlos'][itime,Ibm,:])
            AlldVlos = numpy.ravel(acDict['dVlos'][itime,Ibm,:])
            Allmob = numpy.ravel(acDict['mob'][itime,Ibm,:])
            Allkappa = numpy.ravel(acDict['kappa'][itime,Ibm,:])
            Allchi2 = numpy.ravel(acDict['chi2'][itime,Ibm,:])
            Allfitcode = numpy.ravel(acDict['fitcode'][itime,Ibm,:])




            # setting upper limit at 150 km
            Iht = numpy.where(AllAltitude>htout[-1,-1])[0]
            AllVlos[Iht] = numpy.nan
            AlldVlos[Iht] = numpy.nan

            print('Allchi2', Allchi2)

            # added 10-15-2019 naning out any LOS velocities which are deemed poor fits
            # based on conversation over email with Ashton
            # fitcode between 0-4 are OK, <0 are bad, >5 are bad
            # chi square between 0.1-10 are OK
            # print 'AllVlosB4', AllVlos
            Inan = numpy.where((Allchi2 < self.config['FILTERING']['Chi2Min']) | \
                    (Allchi2 > self.config['FILTERING']['Chi2Max']) | \
                    (Allfitcode < self.config['FILTERING']['FitCodeMin']) | \
                    (Allfitcode > self.config['FILTERING']['FitCodeMax']))

            # Inan = numpy.where((Allfitcode < 1.) | (Allfitcode > 5.))
            AllVlos[Inan] = numpy.nan
            AlldVlos[Inan] = numpy.nan
            # print 'Inan', Inan
            # print 'AC Chi2', Allchi2[Inan]
            # print 'AC Fit codes', Allfitcode[Inan]
            # print 'AllVlos.shape not Nan', AllVlos[~numpy.isnan(AllVlos)].shape
            # print 'AC Nan Shape', Inan[0].shape
            # print 'All Altitude Nan', AllAltitude[Inan]



            '''
            12-02-2019 - after conversation with Ashton for the time being will scale the
            AC dLOS velocities since the Chi 2 appears to be off by approximately 1/sqrt(10)
            Ashton is going to work on fixing this issue and then we can simply set this back.
            '''

            AlldVlos = AlldVlos/numpy.sqrt(self.config['PROCESSING']['dVLOSScaler'])
            # print 'dVlos', AlldVlos
            # # print 'AllVlos - Check', AllVlos
            # # print 'Inan', Inan
            # print '\n \n'

            '''
            added on 18 Feb 2020 - nuin scaler
            if I multiply this to all kappa, it changes the electric field.
            '''

            #print 'nuinScaler', self.nuinScaler
            #Allkappa = Allkappa*self.nuinScaler

            """
            # grabbing corresponding Long Pulse data
            # repeat the data over the interval
            """

            LPUnixTimeMean = numpy.mean(lpDict['UnixTime'], axis=1)
            I = numpy.where((acDict['UnixTime'][itime,0]<=LPUnixTimeMean) \
                        & (acDict['UnixTime'][itime,1]>=LPUnixTimeMean))[0]
            I = numpy.unique(I)

            """
            # assume the beam pattern is the same between LP and AC
            #repeat this over the
            """

            LPBabs = lpDict['Babs'][Ibm,:]
            LPBabs = numpy.ravel(numpy.repeat(LPBabs[numpy.newaxis,:,:],len(I),axis=0))
            LPdec = lpDict['dec'][Ibm,:]
            LPdec = numpy.ravel(numpy.repeat(LPdec[numpy.newaxis,:,:],len(I),axis=0))
            LPdip = lpDict['dip'][Ibm,:]
            LPdip = numpy.ravel(numpy.repeat(LPdip[numpy.newaxis,:,:],len(I),axis=0))
            LPk1 = lpDict['k1'][Ibm,:,:]
            LPk1 = numpy.repeat(LPk1[numpy.newaxis,:,:,:], len(I), axis=0)
            LPk1 = numpy.reshape(LPk1,(Ibm.shape[0]*lpDict['NAltitudes']*len(I),3))
            LPAltitude = lpDict['Altitude'][Ibm,:]
            LPAltitude = numpy.ravel(numpy.repeat(LPAltitude[numpy.newaxis,:,:],len(I),axis=0))/1000.

            # can't seem to do this in one line?
            # give shape mismatch error
            # lpDict['Vlos'][I,Ibm,:] does not work
            LPVlos = lpDict['Vlos'][:,Ibm,:]
            LPVlos = numpy.ravel(LPVlos[I,:,:])
            LPdVlos = lpDict['dVlos'][:,Ibm,:]
            LPdVlos = numpy.ravel(LPdVlos[I,:,:])
            LPmob = lpDict['mob'][:,Ibm,:]
            LPmob = numpy.ravel(LPmob[I,:,:])
            LPkappa = lpDict['kappa'][:,Ibm,:]
            LPkappa = numpy.ravel(LPkappa[I,:,:])
            LPchi2 = lpDict['chi2'][:,Ibm,:]
            LPchi2 = numpy.ravel(LPchi2[I,:,:])
            LPfitcode = lpDict['fitcode'][:,Ibm,:]
            LPfitcode = numpy.ravel(LPfitcode[I,:,:])

            # do the nan filtering for the LP data
            InanLP = numpy.where((LPchi2 < self.config['FILTERING']['Chi2Min']) |\
                    (LPchi2 > self.config['FILTERING']['Chi2Max']) | \
                    (LPfitcode < self.config['FILTERING']['FitCodeMin']) | \
                    (LPfitcode > self.config['FILTERING']['FitCodeMax']))
            #InanLP = numpy.where((LPfitcode < 1.) | (LPfitcode > 5.))
            LPVlos[InanLP] = numpy.nan
            # print 'NanLP', InanLP
            # print 'LPchi2', LPchi2[InanLP]
            # print 'LPfitcode', LPfitcode[InanLP]
            # print 'LPVlos.shape not Nan', LPVlos[~numpy.isnan(LPVlos)].shape
            # print 'NanLP shape', InanLP[0].shape
            #print '\n \n'

            # do altitude filtering
            ACAltFilter = numpy.where((AllAltitude<=htout[-1,-1]) \
                            & (AllAltitude>=htout[0,0]))

            LPAltFilter = numpy.where((LPAltitude>=LPminAlt) & (LPAltitude<=LPmaxAlt))[0]

            # merge the AC code data with the LP data
            AllBabs = numpy.concatenate((AllBabs[ACAltFilter],LPBabs[LPAltFilter]))
            Alldec = numpy.concatenate((Alldec[ACAltFilter],LPdec[LPAltFilter]))
            Alldip = numpy.concatenate((Alldip[ACAltFilter],LPdip[LPAltFilter]))
            Allk = numpy.concatenate((Allk[ACAltFilter],LPk1[LPAltFilter]))
            AllAltitude = numpy.concatenate((AllAltitude[ACAltFilter], LPAltitude[LPAltFilter])) # convert to km
            AllVlos = numpy.concatenate((AllVlos[ACAltFilter], LPVlos[LPAltFilter]))
            AlldVlos = numpy.concatenate((AlldVlos[ACAltFilter], LPdVlos[LPAltFilter]))
            Allmob = numpy.concatenate((Allmob[ACAltFilter], LPmob[LPAltFilter]))
            Allkappa = numpy.concatenate((Allkappa[ACAltFilter], LPkappa[LPAltFilter]))
            # Allk = numpy.concatenate((Allk[ACAltFilter,:],LPk1[LPAltFilter,:]), axis=0)


            # sys.exit()

            # test in the case that kappa >>1, will be collisionless case
            if self.config['PROCESSING']['TestAlgorithm']  == True:
                tempArr = numpy.ones(Allkappa.shape)*1e20
                Allkappa = tempArr.copy()


            # averaging over the declination and the dip angle
            decAltGrid=numpy.zeros(htout.shape[0],dtype='float64')
            dipAltGrid=numpy.zeros(htout.shape[0],dtype='float64')
            BabsAltGrid=numpy.zeros(htout.shape[0],dtype='float64')
            for bb in range(htout.shape[0]):
                I=numpy.where((AllAltitude>=htout[bb,0])&(AllAltitude<=htout[bb,1]))[0]
                decAltGrid[bb]=numpy.nanmean(Alldec[I])
                dipAltGrid[bb]=numpy.nanmean(Alldip[I])
                BabsAltGrid[bb]=numpy.nanmean(AllBabs[I])



            '''
            # invert the winds
            # 09-17-2017 will need to work through this to know what issues are and fix them.
            #htout,numpy.array(test),numpy.array(terr),numpy.array(VlosEst),Iout
            ###
            #NO NANS GO INTO THIS ALGORITHM!!!!!!!!!!!!!!!!!!!!!!!!!

            ###
            '''
            try:
                (htout,EstimatedWinds,\
                CovEstimatedWind,LosVelEstimated,Iout,\
                VlosOut,dVlosOut,AllAltitudeOut) = self.winds.invertwinds(\
                                                        htout,AllAltitude,AllVlos,\
                                                        AlldVlos,Allk,Allmob,Allkappa,\
                                                        decAltGrid,dipAltGrid, \
                                                        abserrFilter = self.config['INVERTWINDS']['AbsoluteErrorFilter'], \
                                                        windCovar = self.config['INVERTWINDS']['Wind_ModelCovariance'], \
                                                        EfieldCovar = self.config['INVERTWINDS']['Efield_ModelCovariance']\
                                                        )
            except:
                print('failed to run invertWinds')
                statusArr[itime] = -1
                continue
                # break
            # diagnostic print statement



            # first 3 elements are the electric field
            # remaining elements are the wind components on the htout grid
            # should be 3+htout.shape[0] elements long
            outDict['Efield'][itime,:] = numpy.transpose(EstimatedWinds[0:3])
            outDict['errEfield'][itime,:] = numpy.sqrt(numpy.diag(EstimatedWinds)[0:3]) # grab the diagonal elements
            tempWind = EstimatedWinds[3:]
            tempErrWind = numpy.sqrt(numpy.diag(CovEstimatedWind)[3:])

            # print 'tempWind.shape', tempWind.shape
            # print 'tempWindErr Shape', tempErrWind.shape
            # print '\n\n'
            # print 'tempWind', tempWind
            # print 'tempErrWind', tempErrWind
            # print '\n\n'
            # print tempWind[0::3],tempErrWind[0::3]
            #
            # print 'VlosEstimated', LosVelEstimated, AllVlos
            # print 'VlosEstimate.shape AllVlos.shape', LosVelEstimated.shape, AllVlos.shape

            # print numpy.ravel(VlosOut), numpy.ravel(LosVelEstimated)
            chi2 = (numpy.ravel(VlosOut)-numpy.ravel(LosVelEstimated))**2/(numpy.ravel(dVlosOut)**2)
            # print numpy.ravel(VlosOut)-numpy.ravel(LosVelEstimated)
            # print dVlosOut, AllAltitudeOut
            # print '\n\n'
            # print 'Chi2', numpy.sum(chi2)/(numpy.ravel(VlosOut).shape[0]-1.)
            # print AllAltitude.shape
            # print '\n\n'

            ## print outDict.keys()
            #NN = numpy.ravel(VlosOut).shape[0]
            #outDict['Vlos'][itime,0:NN] = numpy.ravel(VlosOut)
            #outDict['dVlos'][itime,0:NN] = numpy.ravel(dVlosOut)
            #outDict['VlosEst'][itime,0:NN] = numpy.ravel(LosVelEstimated)

            #print(htout)
            #print(htout.shape[0],htout.shape[1])
            #for iiht in range(htout.shape[0]):
            #    qalt = numpy.where((AllAltitude >= htout[iiht,0]) & (AllAltitude <= htout[iiht,1]) )[0]
            #    print('htout', htout[iiht,0], htout[iiht,1])
            #    print('qalt', qalt)
            #    print(AllVlos[qalt])
            #    NN = numpy.ravel(AllVlos[qalt]).shape[0]
            #    outDict['VlosAltGrid'][itime,iiht,0:NN] = numpy.ravel(AllVlos[qalt])
            #    outDict['dVlosAltGrid'][itime,iiht,0:NN] = numpy.ravel(AlldVlos[qalt])

            #perp east, north, up
            outDict['WindGmag'][itime,:,0] = tempWind[0::3][numpy.newaxis,:,0]
            outDict['WindGmag'][itime,:,1] = tempWind[1::3][numpy.newaxis,:,0]
            outDict['WindGmag'][itime,:,2] = tempWind[2::3][numpy.newaxis,:,0]
            outDict['errWindGmag'][itime,:,0] = tempErrWind[0::3]
            outDict['errWindGmag'][itime,:,1] = tempErrWind[1::3]
            outDict['errWindGmag'][itime,:,2] = tempErrWind[2::3]

            WindGmag = outDict['WindGmag']
            for ii in range(htout.shape[0]):
                outDict['WindGeo'][itime,ii,:] = numpy.squeeze(self.winds.gmag2geo(\
                                                numpy.array([[WindGmag[itime,ii,0]],\
                                                [WindGmag[itime,ii,1]],\
                                                [WindGmag[itime,ii,2]]]),\
                                                numpy.deg2rad(decAltGrid[ii]),\
                                                numpy.deg2rad(dipAltGrid[ii])))
            terrgeo = self.winds.gmag2geo_covar(CovEstimatedWind[3:,3:],\
                                            numpy.deg2rad(decAltGrid),\
                                            numpy.deg2rad(dipAltGrid))
            tempErrWindGeo = numpy.sqrt(numpy.diag(terrgeo))
            outDict['errWindGeo'][itime,:,0] = tempErrWindGeo[0::3]
            outDict['errWindGeo'][itime,:,1] = tempErrWindGeo[1::3]
            outDict['errWindGeo'][itime,:,2] = tempErrWindGeo[2::3]

            # diagnostic print statements
            # print 'Winds 0 component geomag', outDict['WindGmag'][itime,:,0]
            # print 'Winds 0 component geographic', outDict['WindGeo'][itime,:,0]
            # print '-------------------------'
            # print '\n\n'

            # now need to calculate the velocities in the E-region
            # need this to calculate Joule Heating and currents


            # pIn = [self.config['VVELS']['FracErrorOffset'],self.config['VVELS']['FracErrorThreshold'],\
            #         numpy.nan,self.config['VVELS']['AbsoluteErrorThreshold']] # parameters in
                    # nan position is unused, but I don't want to change the algorithm
            # PLAT,AllVlos,AlldVlos,Allk,AllPlat,AllPlong,Allht,htmin=150.0*1000,htmax=400.0*1000,\
        	    # covar=[1000.*1000.,1000.*1000.,5.*5.],p=[200.0,0.5,2000.0,100.0])
            (alt_out1,tempVest,tempdVest,\
            tempdVestAll,Foo,statusOut) = self.winds.compute_velvec2(htout,\
                                                AllVlos,AlldVlos,Allk,\
                                                AllAltitude, [],AllAltitude,\
                                                htmin = self.config['VVELS']['Velocity_MinAltitude'],\
                                                htmax = self.config['VVELS']['Velocity_MaxAltitude'],\
                                                covar = self.config['VVELS']['Velocity_ModelCovariance'],\
                                                FracErrorOffset = self.config['VVELS']['FracErrorOffset'], \
                                                FracErrorThreshold =  self.config['VVELS']['FracErrorThreshold'], \
                                                AbsoluteErrorThreshold = self.config['VVELS']['AbsoluteErrorThreshold'])

            if statusOut == False:
                statusArr[itime] = -2

            ## print 'tmpVest.shape', tempVest.shape
            #outDict['VestGmag'][itime,:,:] = tempVest
            #outDict['errVestGmag'][itime,:,:] = tempdVest

            ## note this is different from what Mike N Did
            ## you don't which beams were chosen, so you need to loop through the chosen
            ## beams and get their index correct

            #tmphtout = numpy.array([[200.,400.]])
            #(alt_out1,tempVest,tempdVest,\
            #tempdVestAll,Foo,statusOut) = self.winds.compute_velvec2(tmphtout,\
            #                                    AllVlos,AlldVlos,Allk,\
            #                                    AllAltitude, [],AllAltitude,\
            #                                    htmin = self.config['VVELS']['Velocity_MinAltitude'],\
            #                                    htmax = self.config['VVELS']['Velocity_MaxAltitude'],\
            #                                    covar = self.config['VVELS']['Velocity_ModelCovariance'],\
            #                                    FracErrorOffset = self.config['VVELS']['FracErrorOffset'], \
            #                                    FracErrorThreshold =  self.config['VVELS']['FracErrorThreshold'], \
            #                                    AbsoluteErrorThreshold = self.config['VVELS']['AbsoluteErrorThreshold'])
            ##
            #if statusOut == False:
            #    statusArr[itime] = -3

            #outDict['VestGmag_300km'][itime,:] = tempVest
            #outDict['errVestGmag_300km'][itime,:] = tempdVest
            #dotprod = outDict['VestGmag'][itime,:,0]*outDict['VestGmag_300km'][itime,0]+\
            #        outDict['VestGmag'][itime,:,1]*outDict['VestGmag_300km'][itime,1]+\
            #        outDict['VestGmag'][itime,:,2]*outDict['VestGmag_300km'][itime,2]
            #normVest = numpy.sqrt(outDict['VestGmag'][itime,:,0]**2 + outDict['VestGmag'][itime,:,1]**2 + outDict['VestGmag'][itime,:,2]**2)
            #normV300 = numpy.sqrt(outDict['VestGmag_300km'][itime,0]**2 + outDict['VestGmag_300km'][itime,1]**2 + outDict['VestGmag_300km'][itime,2]**2)
            #outDict['Angle'][itime,:] = numpy.rad2deg(numpy.arccos(dotprod/(normVest*normV300)))


            #"""
            #Coriolis and Centrifugal terms
            #"""
            ## calculate the colatitude
            #ThetaRadian = numpy.deg2rad(90.-acDict['SiteLatitude'])
            ## radius array
            #rArr = (6371.+htoutm)*1000. # meters
            ## outDict['WindGeo'][itime,ii,:]
            ## derived this expression and it is in Fuller-Rowell and Rees, 1984
            #AcentrifugalN = -(outDict['WindGeo'][itime,:,0]**2.)*numpy.cos(ThetaRadian)/(rArr*numpy.sin(ThetaRadian))
            #AcentrifugalU = -(outDict['WindGeo'][itime,:,0]**2.)/rArr

            #outDict['Centrifugal'][itime,:,1] = AcentrifugalN
            #outDict['Centrifugal'][itime,:,2] = AcentrifugalU

            ## angular velocity of earth's rotation
            #OmegaE = 7.2921150e-5 #rads/s
            #AcoriolisE = outDict['WindGeo'][itime,:,1]*OmegaE*numpy.cos(ThetaRadian)
            #AcoriolisN = -outDict['WindGeo'][itime,:,0]*OmegaE*numpy.cos(ThetaRadian)
            #outDict['Coriolis'][itime,:,0] = AcoriolisE
            #outDict['Coriolis'][itime,:,1] = AcoriolisN

            #"""
            #The Joule Heating Calculations
            #"""
            #for ibeam in Ibm:
            #    #htoutm = height out mean
            #    Ve_onNeGrid = interp1d(htoutm,outDict['VestGmag'][itime,:,0],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
            #    Vn_onNeGrid = interp1d(htoutm,outDict['VestGmag'][itime,:,1],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
            #    Ue_onNeGrid = interp1d(htoutm,outDict['WindGmag'][itime,:,0],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
            #    Un_onNeGrid = interp1d(htoutm,outDict['WindGmag'][itime,:,1],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
            #    tempB = interp1d(htoutm,BabsAltGrid,bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)

                ## print 'Vn_onNeGrid',  Vn_onNeGrid, Ve_onNeGrid, outDict['Efield'][itime,1]/tempB, outDict['Efield'][itime,0]/tempB
                ## for ialt in range(acDict['Altitude'].shape[-1]):
                ##     print 'Ue_onNeGrid', Ue_onNeGrid[ialt],acDict['Altitude'][ibeam,ialt]/1000.
                ## sys.exit()

                ## fill the variables
                ## see equation Thayer 1998, JGR, equation 8
                ##east = x, north = y
                #outDict['CurrentGmag'][itime,ibeam,:,0] = elementary_charge*acDict['Ne'][itime,ibeam,:]*(Ve_onNeGrid + outDict['Efield'][itime,1]/tempB)
                #outDict['CurrentGmag'][itime,ibeam,:,1] = elementary_charge*acDict['Ne'][itime,ibeam,:]*(Vn_onNeGrid - outDict['Efield'][itime,0]/tempB)

                ## print outDict['CurrentGmag'][:,ibeam]

                ## joule heating, electric field only
                #outDict['JouleHeatingE'][itime,ibeam,:] = acDict['PedersenConductivity'][itime,ibeam,:]*(outDict['Efield'][itime,:]**2).sum()
                #indx = numpy.where(numpy.isfinite(outDict['JouleHeatingE'][itime,ibeam,:]) &\
                #                    numpy.isfinite(acDict['Altitude'][ibeam,:]) )[0]
                #outDict['IntegratedJouleHeatingE'][itime,ibeam] = trapz(outDict['JouleHeatingE'][itime,ibeam,indx],acDict['Altitude'][ibeam,indx])

                ## Joule heating with neutral wind
                ## see equation Thayer 1998, JGR, equation 1
                ## minor sign error fixed on Ue on 11/14/2018
                ## to be consistent with equation 1
                #outDict['JouleHeatingTotal'][itime,ibeam,:] = acDict['PedersenConductivity'][itime,ibeam,:]*\
                #                                            ((outDict['Efield'][itime,0] - Un_onNeGrid*tempB)**2 + \
                #                                            (outDict['Efield'][itime,1] + Ue_onNeGrid*tempB)**2)
                #indx = numpy.where(numpy.isfinite(outDict['JouleHeatingTotal'][itime,ibeam,:]) &\
                #                    numpy.isfinite(acDict['Altitude'][ibeam,:]) )[0]
                #outDict['IntegratedJouleHeatingTotal'][itime,ibeam] = trapz(outDict['JouleHeatingTotal'][itime,ibeam,indx],acDict['Altitude'][ibeam,indx])

                #### add the mechanical heating terms
                ## added on 07 01 2020
                ## see notes.
                ## -jNorth*Ueast*B
                #UnJCrossB_East = -outDict['CurrentGmag'][itime,ibeam,:,1]*Ue_onNeGrid*tempB

                ## jEast*Unorth*tempB
                #UnJCrossB_North = outDict['CurrentGmag'][itime,ibeam,:,0]*Un_onNeGrid*tempB
                #tmpTotalMechanicalHeating = UnJCrossB_East+UnJCrossB_North

                #outDict['JouleHeatingMechanical'][itime,ibeam,:] = tmpTotalMechanicalHeating
                #indx = numpy.where(numpy.isfinite(outDict['JouleHeatingMechanical'][itime,ibeam,:]) &\
                #                    numpy.isfinite(acDict['Altitude'][ibeam,:]) )[0]
                #outDict['IntegratedJouleHeatingMechanical'][itime,ibeam] = trapz(outDict['JouleHeatingMechanical'][itime,ibeam,indx], acDict['Altitude'][ibeam,indx])


                #"""
                #07 26 2021
                #calculate the other thayer terms
                #"""
                ## print 'winds gmag shape', outDict['WindGmag'][itime,:,0].shape
                ## print 'current gmag shape outDict', outDict['CurrentGmag'].shape
                ## print 'tmpB,', tempB.shape
                ## print 'tmpkappa,', acDict['kappa'][itime,ibeam,:]
                ## print 'BabsAltGrid', BabsAltGrid.shape
                ## print 'Babs acDict', acDict['Babs'][ibeam,:]
                ## print 'Ne', acDict['Ne'][itime,ibeam,:]
                ## print 'tmpkappa', acDict['kappa'][0,0,:]
                ## BabsFullGrid = numpy.tile(acDict['Babs'], acDict['Ne'].shape[0]).reshape([acDict['Ne'].shape[0], acDict['Ne'].shape[1], acDict['Ne'].shape[2]])
                ## print 'BabsFullGrid shape', BabsFullGrid.shape
                ## equation 10 in Thayer 1998
                ## print outDict['CurrentGmag'][0,-1,:]

                ## print 'j2.shape', j2.shape
                ## print 'j2', j2

                #"""
                #Thayer 2000 JGR Table 1
                #"""
                #j2 = outDict['CurrentGmag'][itime,ibeam,:,0]**2+outDict['CurrentGmag'][itime,ibeam,:,1]**2
                #outDict['JouleHeatingTotalThayer'][itime,ibeam,:] = j2*acDict['Babs'][ibeam,:]*acDict['kappa'][itime,ibeam,:]/(elementary_charge*acDict['Ne'][itime,ibeam,:])

                ## j dot E from table
                #tmpQem = outDict['CurrentGmag'][itime,ibeam,:,0]*outDict['Efield'][itime,0] + outDict['CurrentGmag'][itime,ibeam,:,1]*outDict['Efield'][itime,1]
                #outDict['EMTransferRateThayer'][itime,ibeam,:] = tmpQem

                ## q-qJ from Thayer 2000 Table 1
                #outDict['JouleHeatingMechanicalThayer'][itime,ibeam,:] = tmpQem - outDict['JouleHeatingTotalThayer'][itime,ibeam,:]
                ## print JouleHeatingTotalThayer.shape
                ## print outDict['JouleHeatingTotal'].shape
                ## print outDict['JouleHeatingTotalThayer'].shape


            # print "Joule heating difference"
            # for itime in range(47):
            #     print numpy.nanmax((outDict['JouleHeatingTotal'][itime,:] - outDict['JouleHeatingTotalThayer'][itime,:])/outDict['JouleHeatingTotal'][itime,:])
            # sys.exit()




            """
            Anything along the vertical beam
            """
            # now grab the vertical beam
            # these are on the ISR experiment resolution - need to interpolate onto standard grid
            # edited on 28 January 2020
            if len(qvert) != 0:

                tempVertAlt = numpy.ravel(acDict['Altitude'][qvert,:]/1000.)
                # print tempVertAlt.shape,numpy.ravel(outDict['JouleHeatingE'][itime,qvert,:])
                #tmpJoule = interp1d(tempVertAlt,numpy.ravel(outDict['JouleHeatingE'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamJouleHeatingE'][itime,:] = tmpJoule
                #tmpJouleTot = interp1d(tempVertAlt,numpy.ravel(outDict['JouleHeatingTotal'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamJouleHeatingTotal'][itime,:] = tmpJouleTot
                #tmpJouleMech = interp1d(tempVertAlt,numpy.ravel(outDict['JouleHeatingMechanical'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamJouleHeatingMechanical'][itime,:] = tmpJouleMech

                # added on 05/03/2022 - try to get full altitude




                # adding in Thayer's stuff.
                """
                Adding in Thayer vertical beam
                """
                #tmpJouleTot_Thayer = interp1d(tempVertAlt,numpy.ravel(outDict['JouleHeatingTotalThayer'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamJouleHeatingTotal_Thayer'][itime,:] = tmpJouleTot_Thayer

                #tmpEMRate_Thayer = interp1d(tempVertAlt,numpy.ravel(outDict['EMTransferRateThayer'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamEMTranfer_Thayer'][itime,:] = tmpEMRate_Thayer

                #tmpMechHeating_Thayer = interp1d(tempVertAlt,numpy.ravel(outDict['JouleHeatingMechanicalThayer'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamJouleHeatingMechanical_Thayer'][itime,:] = tmpMechHeating_Thayer



                # adding in the Lorentz forcing
                #interpolate the Hall and Pedersen conductivity onto the grid
                #BabsAltGrid is on the altitude grid and probbly good enough
                #tmpHallCond = interp1d(tempVertAlt,numpy.ravel(outDict['HallConductivity'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpPedCond = interp1d(tempVertAlt,numpy.ravel(outDict['PedersenConductivity'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['VerticalBeamHallConductivity'][itime,:] = tmpHallCond
                #outDict['VerticalBeamPedersenConductivity'][itime,:] = tmpPedCond

                ## 05/03/2022
                #FullAltitude = numpy.arange(90., 255.,5.)
                #tmpPedCondFullAlt = interp1d(tempVertAlt,numpy.ravel(outDict['PedersenConductivity'][itime,qvert,:]),bounds_error=0,fill_value=numpy.nan)(FullAltitude)
                #E2 = outDict['Efield'][itime,0]**2 + outDict['Efield'][itime,1]**2
                #tmpJHFull = tmpPedCondFullAlt*E2
                #IntegratedJHFull = trapz(tmpJHFull,FullAltitude*1000.)
                #IntegratedPedConductFull = trapz(tmpPedCondFullAlt,FullAltitude*1000.)

                #outDict['AltitudeFull'] = FullAltitude
                #outDict['PedersenConductivtityFull'][itime,:] = tmpPedCondFullAlt
                #outDict['JouleHeatingEFull'][itime,:] = tmpJHFull
                #outDict['PedersenConductanceFull'][itime] = IntegratedPedConductFull
                ##outDict['IntegratedJouleHeatingEFull'][itime,:] = IntegratedJHFull
                #outDict['IntegratedJouleHeatingEFull'][itime] = IntegratedJHFull


                ##tmpZonalFlow = -MeridE/0.495595*1e4 #tesla
                ##tmpMeridFlow = ZonalE/0.495595*1e4 # convert to Tesla
                #VplasmaE = -outDict['Efield'][itime,1]/0.495595*1e4
                #VplasmaN = outDict['Efield'][itime,0]/0.495595*1e4
                #nMass = interp1d(tempVertAlt,numpy.ravel(acDict['nMass'][itime,qvert,:]),bounds_error=0)(htoutm)

                #print('tmpHallCond', tmpHallCond)
                #print('BabsAltGrid', BabsAltGrid)
                #print('nmass', nMass)

                #alphaH = tmpHallCond*BabsAltGrid*BabsAltGrid/nMass
                #alphaP = tmpPedCond*BabsAltGrid*BabsAltGrid/nMass

                #print('alphaP', alphaP)
                #print('alphaH', alphaH)

                ##outDict['HallDrag'][itime,:] = alphaH
                ##outDict['PedersenDrag'][itime,:] = alphaP

                ## Larsen Waltersheid 1995 eq 5
                #FLorentzE = alphaH*(VplasmaN-outDict['WindGmag'][itime,:,1]) + alphaP*(VplasmaE-outDict['WindGmag'][itime,:,0])
                #FLorentzN = alphaP*(VplasmaN-outDict['WindGmag'][itime,:,1]) - alphaH*(VplasmaE-outDict['WindGmag'][itime,:,0])
                ##outDict['Lorentz'][itime,:,0] = FLorentzE
                ##outDict['Lorentz'][itime,:,1] = FLorentzN


                #print('FLorentzN', FLorentzN)


                #tmpkappa = interp1d(tempVertAlt,numpy.ravel(acDict['kappa'][itime,qvert,:]),bounds_error=0)(htoutm)
                #outDict['Kappa'][itime,:] = tmpkappa

                #tmpnuin = interp1d(tempVertAlt,numpy.ravel(acDict['nuin'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpnuinBrekke = interp1d(tempVertAlt,numpy.ravel(acDict['nuinBrekke'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpTi = interp1d(tempVertAlt,numpy.ravel(acDict['Ti'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpdTi = interp1d(tempVertAlt,numpy.ravel(acDict['dTi'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpTn = interp1d(tempVertAlt,numpy.ravel(acDict['Tn'][itime,qvert,:]),bounds_error=0)(htoutm)
                #print('Kappa', tmpkappa)
                #print('\n\n')

                ##outDict['Vertical_nuin'][itime,:] = tmpnuin
                ##outDict['Vertical_nuin_Brekke'][itime,:] = tmpnuinBrekke
                ##outDict['VerticalTi'][itime,:] = tmpTi
                ##outDict['VerticaldTi'][itime,:] = tmpdTi
                ##outDict['VerticalTn'][itime,:] = tmpTn

                #tmpNe = interp1d(tempVertAlt,numpy.ravel(acDict['Ne'][itime,qvert,:]),bounds_error=0)(htoutm)
                #tmpdNe = interp1d(tempVertAlt,numpy.ravel(acDict['dNe'][itime,qvert,:]),bounds_error=0)(htoutm)
                ##outDict['VerticalBeamNe'][itime,:] = tmpNe
                ##outDict['errVerticalBeamNe'][itime,:] = tmpdNe


            ## average over altitude grid
            #tmpAltitude = numpy.ravel(acDict['Altitude'])/1000.
            #tmpJouleE = numpy.ravel(outDict['JouleHeatingE'][itime,:,:])
            #tmpJouleTotal = numpy.ravel(outDict['JouleHeatingTotal'][itime,:,:])
            #tmpJouleMech = numpy.ravel(outDict['JouleHeatingMechanical'][itime,:,:])
            #tmpHallCond = numpy.ravel(outDict['HallConductivity'][itime,:,:])
            #tmpPedCond = numpy.ravel(outDict['PedersenConductivity'][itime,:,:])
            #outDict['AltitudeJHBeam'] = acDict['Altitude']

            ### add in the Thayer information
            ### outDict['EMTransferRateThayer'], outDict['JouleHeatingMechanicalThayer']
            ##tmpJouleTotal_Thayer = numpy.ravel(outDict['JouleHeatingTotalThayer'][itime,:,:])
            ##tmpEMTransfer_Thayer = numpy.ravel(outDict['EMTransferRateThayer'][itime,:,:])
            ##tmpJouleMech_Thayer = numpy.ravel(outDict['JouleHeatingMechanicalThayer'][itime,:,:])

            ## adding in SNR and Ne - unfitted
            #tmpAltitudeRaw = numpy.ravel(acDict['AltitudeRaw'])/1000.
            #tmpSNRRaw = numpy.ravel(acDict['SNRRaw'][itime,:,:])
            #tmpNeFitted = numpy.ravel(acDict['Ne'][itime,:,:]) # fitted
            #tmpNeRaw = numpy.ravel(acDict['NeRaw'][itime,:,:])

            #tmpScaleHeight = numpy.ravel(acDict['ScaleHeight'][itime,:,:])

            #for bb in range(htout.shape[0]):
            #    # putting everything onto the same altitude grid as the winds
            #    I = numpy.where((tmpAltitude >= htout[bb,0]) & (tmpAltitude <= htout[bb,1]))[0]
            #    IRaw = numpy.where((tmpAltitudeRaw >= htout[bb,0]) & (tmpAltitudeRaw <= htout[bb,1]))[0]
            #    Iend = I.shape[0]
            #    IendRaw = IRaw.shape[0]
            #    #outDict['MedianJouleHeatingE'][itime,bb] = numpy.nanmedian(tmpJouleE[I])
            #    #outDict['MedianJouleHeatingTotal'][itime,bb] = numpy.nanmedian(tmpJouleTotal[I])
            #    #outDict['MedianJouleHeatingMechanical'][itime,bb] = numpy.nanmedian(tmpJouleMech[I])
            #    #outDict['MedianHallConductivity'][itime,bb] = numpy.nanmedian(tmpHallCond[I])
            #    #outDict['MedianPedersenConductivity'][itime,bb] = numpy.nanmedian(tmpPedCond[I])

            #    ## add in thayer stuff.
            #    #outDict['MedianJouleHeatingTotal_Thayer'][itime,bb] = numpy.nanmedian(tmpJouleTotal_Thayer[I])
            #    #outDict['MedianJouleHeatingMechanical_Thayer'][itime,bb] = numpy.nanmedian(tmpJouleMech_Thayer[I])
            #    #outDict['MedianEMTransfer_Thayer'][itime,bb] = numpy.nanmedian(tmpEMTransfer_Thayer[I])


            #    #outDict['MeanJouleHeatingE'][itime,bb] = numpy.nanmean(tmpJouleE[I])
            #    #outDict['MeanJouleHeatingTotal'][itime,bb] = numpy.nanmean(tmpJouleTotal[I])
            #    #outDict['MeanJouleHeatingMechanical'][itime,bb] = numpy.nanmean(tmpJouleMech[I])
            #    #outDict['MeanHallConductivity'][itime,bb] = numpy.nanmean(tmpHallCond[I])
            #    #outDict['MeanPedersenConductivity'][itime,bb] = numpy.nanmean(tmpPedCond[I])

            #    #outDict['MeanJouleHeatingTotal_Thayer'][itime,bb] = numpy.nanmean(tmpJouleTotal_Thayer[I])
            #    #outDict['MeanJouleHeatingMechanical_Thayer'][itime,bb] = numpy.nanmean(tmpJouleMech_Thayer[I])
            #    #outDict['MeanEMTransfer_Thayer'][itime,bb] = numpy.nanmean(tmpEMTransfer_Thayer[I])


            #    #outDict['MeanNeFitted'][itime,bb] = numpy.nanmean(tmpNeFitted[I])
            #    #outDict['MedianNeFitted'][itime,bb] = numpy.nanmedian(tmpNeFitted[I])
            #    #outDict['MeanSNR'][itime,bb] = numpy.nanmean(tmpSNRRaw[IRaw])
            #    #outDict['MedianSNR'][itime,bb] = numpy.nanmedian(tmpSNRRaw[IRaw])
            #    #outDict['MeanNeRaw'][itime,bb] = numpy.nanmean(tmpNeRaw[IRaw])
            #    #outDict['MedianNeRaw'][itime,bb] = numpy.nanmedian(tmpNeRaw[IRaw])

            #    #outDict['ScaleHeight'][itime,bb] = numpy.nanmean(tmpScaleHeight[I])

            #    #outDict['NeFittedRaw'][itime,bb,0:Iend] = tmpNeFitted[I]
            #    #outDict['SNRRaw'][itime,bb,0:IendRaw] = tmpSNRRaw[IRaw]



            ## now integrate to get
            #tmpAltitude = numpy.nanmean(htout, axis=1)*1000. # convert to m
            #outDict['IntegratedMedianJouleHeatingE'][itime] = trapz(outDict['MedianJouleHeatingE'][itime,:],tmpAltitude)
            #outDict['IntegratedMedianJouleHeatingTotal'][itime] = trapz(outDict['MedianJouleHeatingTotal'][itime,:],tmpAltitude)
            #outDict['IntegratedMedianJouleHeatingMechanical'][itime] = trapz(outDict['MedianJouleHeatingMechanical'][itime,:], tmpAltitude)
            #outDict['MedianHallConductance'][itime] = trapz(outDict['MedianHallConductivity'][itime,:], tmpAltitude)
            #outDict['MedianPedersenConductance'][itime] = trapz(outDict['MedianPedersenConductivity'][itime,:], tmpAltitude)

            #outDict['IntegratedMedianJouleHeatingTotal_Thayer'][itime] = trapz(outDict['MedianJouleHeatingTotal_Thayer'][itime,:],tmpAltitude)
            #outDict['IntegratedMedianJouleHeatingMechanical_Thayer'][itime] = trapz(outDict['MedianJouleHeatingMechanical_Thayer'][itime,:],tmpAltitude)
            #outDict['IntegratedMedianEMTransfer_Thayer'][itime] = trapz(outDict['MedianEMTransfer_Thayer'][itime,:],tmpAltitude)

            #outDict['IntegratedMeanJouleHeatingE'][itime] = trapz(outDict['MeanJouleHeatingE'][itime,:],tmpAltitude)
            #outDict['IntegratedMeanJouleHeatingTotal'][itime] = trapz(outDict['MeanJouleHeatingTotal'][itime,:],tmpAltitude)
            #outDict['IntegratedMeanJouleHeatingMechanical'][itime] = trapz(outDict['MeanJouleHeatingMechanical'][itime,:], tmpAltitude)
            #outDict['MeanHallConductance'][itime] = trapz(outDict['MeanHallConductivity'][itime,:], tmpAltitude)
            #outDict['MeanPedersenConductance'][itime] = trapz(outDict['MeanPedersenConductivity'][itime,:], tmpAltitude)

            #outDict['IntegratedMeanJouleHeatingTotal_Thayer'][itime] = trapz(outDict['MeanJouleHeatingTotal_Thayer'][itime,:],tmpAltitude)
            #outDict['IntegratedMeanJouleHeatingMechanical_Thayer'][itime] = trapz(outDict['MeanJouleHeatingMechanical_Thayer'][itime,:],tmpAltitude)
            #outDict['IntegratedMeanEMTransfer_Thayer'][itime] = trapz(outDict['MeanEMTransfer_Thayer'][itime,:],tmpAltitude)

            outDict['nuinScaler'][itime] = acDict['nuinScaler'][itime]
                # fill status
            if statusArr[itime] >= 0:
                statusArr[itime] = 1 # complete status

        # for itime in range(outDict['JouleHeatingTotal'].shape[0]):
        #     for ibeam in Ibm:
        #         for ialt in range(outDict['JouleHeatingTotal'].shape[-1]):
        #             print acDict['Altitude'][ibeam,ialt]/1000, \
        #                     outDict['JouleHeatingTotal'][itime,ibeam,ialt], \
        #                     outDict['JouleHeatingTotalThayer'][itime,ibeam,ialt], \
        #                     (outDict['JouleHeatingTotal'][itime,ibeam,ialt]-outDict['JouleHeatingTotalThayer'][itime,ibeam,ialt])/outDict['JouleHeatingTotal'][itime,ibeam,ialt]

        # sys.exit()

        # filling the output geophysical array
        UnixTimeMean = numpy.mean(acDict['UnixTime'], axis=1)

        if self.config['GEOPHYSICALPARAMETERS']['ByPass'] == True:
            print('bypassing geophys params')
        else:
            KP_AP_Dict = self.geoparams.GetGeophysicalParameters(UnixTimeMean)
            for ikey in KP_AP_Dict.keys():
                outDict[ikey] = KP_AP_Dict[ikey]

                magDict = self.geoparams.ProcessGroundMag(UnixTimeMean[0],UnixTimeMean[-1])
            for ikey in magDict.keys():
                outDict[ikey] = magDict[ikey]

        # filling status and version info into the H5 file
        outDict['Status'] = statusArr
        outDict['Version'] = VersionNumber
        now = datetime.datetime.now()
        outDict['FileCreationDateTime'] = now.strftime('%Y/%m/%d %H:%M:%S')
        print(fname_ac.split('/'))
        #unsure if this is correct
        outDict['ExperimentName'] = fname_ac.split('/')[-1]
        outDict['ExperimentDirectory'] = fname_ac.split('/')[0]
        outDict['LongPulseFile'] = fname_lp
        outDict['AlternatingCodeFile'] = fname_ac
        outDict['ConfigFile'] = self.configFile

        self.create_plots(outDict)

        nuInNumber = float(self.nuinScaler)
        #oname_str = (FnameOut+'_winds_'+'nuin%0.1f_'+VersionNumber+'.h5')%nuInNumber
        OutLocation = os.path.join(FullPathOut, FnameOut)
        print('\n OutLocation', OutLocation)
        self.io.SaveOutDict(OutLocation,acDict,outDict)
            # # print 'errWindGeo', outDict['errWindGeo']
            # print 'Ibm', Ibm
        print('\nCOMPLETE\n')
        os.remove(lock_location)
        return outDict, acDict


    def create_plots(self, outDict):

        self.plotprefix='temp_'
        #os.makedirs(os.path.abspath(self.plotsavedir),exist_ok=True)
        os.makedirs('temp_plots',exist_ok=True)

        # break up arrays into chunks of time no bigger than 24 hours
        chunks_to_plot = list()

        num_times = len(outDict['UnixTime'])
        start_ind = 0
        start_time = outDict['UnixTime'][0,0]
        for i,time_pair in enumerate(outDict['UnixTime']):
            temp_start_time, temp_end_time = time_pair
            time_diff = temp_end_time - start_time
            # Add chunk if over 24 hours elapsed
            if (time_diff >= 24*3600):
                chunks_to_plot.append([start_ind,i])
                start_ind = i
                start_time = temp_start_time
        chunks_to_plot.append([start_ind, num_times])

        num_chunks = len(chunks_to_plot)
        for t, [start_ind,end_ind] in enumerate(chunks_to_plot):
            # if only 1 day worth of data, set t=None so we don't have a
            # 'byDay' in the plot file names
            if (num_chunks == 1):
                #vcom_fname = '{}vvelsnw_vel_comp.png'.format(self.plotprefix)
                wcom_fname = '{}vvelsnw_winds_comp.png'.format(self.plotprefix)
                #vmag_fname = '{}vvelsnw_vel_mag.png'.format(self.plotprefix)
                wmag_fname = '{}vvelsnw_winds_mag.png'.format(self.plotprefix)
            else:
                #vcom_fname = '{}vvelsnw_vel_comp_{}.png'.format(self.plotprefix, t)
                wcom_fname = '{}vvelsnw_winds_comp_{}.png'.format(self.plotprefix, t)
                #vmag_fname = '{}vvelsnw_vel_mag_{}.png'.format(self.plotprefix, t)
                wmag_fname = '{}vvelsnw_winds_mag_{}.png'.format(self.plotprefix, t)

            # make vector plots
            times = outDict['UnixTime'][start_ind:end_ind,:]

            #vels = self.Velocity[start_ind:end_ind,:]
            #covvels = self.VelocityCovariance[start_ind:end_ind,:]
            winds = outDict['WindGeo'][start_ind:end_ind,:]*1000.
            #covwinds = self.ElectricFieldCovariance[start_ind:end_ind,:]*1000.*1000.
            covwinds = numpy.ones(outDict['WindGeo'].shape+(3,))

            summary_plots.plot_components(times, outDict['Altitude'], winds, covwinds,
                            titles=['UE (m/s)','UN (m/s)','UU (m/s)'],
                            ylabel='Alt', clim=[[-150.,150.], [0.,35.]],
                            cmap=['coolwarm', 'turbo'],
                            filename=os.path.join('temp_plots',wcom_fname), scale_factors=[1,1,10])

#            summary_plots.plot_components(times, self.bin_mlat, efs, covefs,
#                            titles=['Ed1 (mV/m)','Ed2 (mV/m)','Ed3 (mV/m)'],
#                            ylabel='Apex MLAT', clim=[[-75., 75.], [0., 15.]],
#                            cmap=['coolwarm', 'turbo'],
#                            filename=os.path.join(self.plotsavedir,ecom_fname), scale_factors=[1,1,10])



            ## make magnitude plots
            ## find index of altitude bin that is closest to alt
            #i = np.argmin(np.abs(self.bin_galt[:,0]-alt))
            #vmag = self.Vgd_mag[start_ind:end_ind,i,:]
            #dvmag = self.Vgd_mag_err[start_ind:end_ind,i,:]
            #vdir = self.Vgd_dir[start_ind:end_ind,i,:]
            #dvdir = self.Vgd_dir_err[start_ind:end_ind,i,:]
            #emag = self.Egd_mag[start_ind:end_ind,i,:]*1000.
            #demag = self.Egd_mag_err[start_ind:end_ind,i,:]*1000.
            #edir = self.Egd_dir[start_ind:end_ind,i,:]
            #dedir = self.Egd_dir_err[start_ind:end_ind,i,:]
            #chi2 = self.ChiSquared[start_ind:end_ind,:]

            #titles = ['V mag. (m/s)', 'V mag. err. (m/s)', 'V dir. (deg)', 'V dir. err. (deg)', '']
            #clim = [[0.,1500.],[0., 350.],[-180., 180.],[0., 35.]]
            #cmap = ['viridis', 'turbo', 'twilight', 'turbo']

            #summary_plots.plot_magnitude(times, self.bin_mlat, vmag, dvmag, vdir, dvdir, chi2,
            #                err_thres=100., mag_thres=100., titles=titles,
            #                ylabel='Apex MLAT', clim=clim, cmap=cmap,
            #                filename=os.path.join(self.plotsavedir,vmag_fname))

            #titles = ['E mag. (mV/m)', 'E mag err (mV/m)', 'E dir (deg)', 'E dir err (deg)', '']
            #clim = [[0.,75.],[0., 15.],[-180., 180.],[0., 35.]]
            #cmap = ['viridis', 'turbo', 'twilight', 'turbo']

            #summary_plots.plot_magnitude(times, self.bin_mlat, emag, demag, edir, dedir, chi2,
            #                err_thres=5., mag_thres=5., titles=titles,
            #                ylabel='Apex MLAT', clim=clim, cmap=cmap,
            #                filename=os.path.join(self.plotsavedir,emag_fname))





def main():

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('config_file', nargs='+', help='A configuration file.')

    args = vars(parser.parse_args())

    #acFile = None
    #lpFile = None
    #ConfigFile = None

    #if sys.argv[1]:
    #    if os.path.isfile(sys.argv[1]):
    #        acFile = sys.argv[1]
    #    else:
    #        raise Exception("Alternating Code file does not exist")
    #else:
    #    raise Exception("No Alternating code file given as argument")

    #if sys.argv[2]:
    #    if os.path.isfile(sys.argv[2]):
    #        lpFile = sys.argv[2]
    #    else:
    #        raise Exception("Long Pulse file does not exist")
    #else:
    #    raise Exception("No Long Pulse file given as argument")

    #if sys.argv[3]:
    #    if os.path.isfile(sys.argv[3]):
    #        ConfigFile = sys.argv[3]
    #    else:
    #        raise Exception("Config File does not exist")
    #else:
    #    raise Exception("No config file given as argument")


    #if ConfigFile and acFile and lpFile:
    #    mainClass = ProcessEregionNeutralWinds(ConfigFile)
    #    out = mainClass.run(acFile,lpFile)

    neutwind = ProcessEregionNeutralWinds(args['config_file'])
    neutwind.run()


if __name__ == "__main__":
    main()


#main(fname_ac,fname_lp,oname)


    # outDict[]

    # Efield[itm,:]=scipy.transpose(test[0:3])
    # errEfield[itm,:]=scipy.sqrt(scipy.diag(terr)[0:3])
    # twind=test[3:]
    # errWind=scipy.sqrt(scipy.diag(terr)[3:])
