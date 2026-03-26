"""
IMPORTANT THIS IS SOME DEPRECATED VERSION OF THE CODE DO NOT USE THIS CODEself.

The main function is ProcessEregionNeutralWinds.py
"""

import IOEregionWinds.IOEregionWinds as IOEregionWinds
import CalculateEregionWinds.CalculateEregionWinds as CalculateEregionWinds
import numpy
from scipy.interpolate import interp1d
from scipy.constants import elementary_charge
from tools.loggerinit.LoggerInit import *
from tools.configreader.ConfigReader import *
import os

#instantiation of classes
io = IOEregionWinds.IOEregionWinds()
winds = CalculateEregionWinds.CalculateEregionWinds()


class ProcessEregionNeutralWinds:

    def __init__(self, configFile = None):
        # check if config file exists
        if configFile:
            if os.path.isfile(configFile):
                self.config = self.ConfigReader.read(configFile)
            else:
                raise Exception('Error: Check config file path')

            # parse config file at this point

        else:
            #if no config file specified use standard inputs
            # parse the config file
            # set the neutral wind altitude grid
            self.NW_MinAltitude = 80.0
            self.NW_MaxAltitude = 150.0
            self.NW_DeltaAltitude = 5.0

            #set long pulse altitude range
            self.LP_MinAltitude = 200.0
            self.LP_MaxAltitude = 400.0
        return

    def _config_parser(self):
        """
        This is a private function designed to parse the config file
        """

        return

    def main(self,fname_ac, fname_lp,oname):
        print 'main'

        # make sure that the files exist
        if os.path.isfile(fname_ac):
            acDict= io.MakeInputDictionary(fname_ac)
        else:
            raise Exception('Error: Alternating Code Input File does not exist')

        if os.path.isfile(fname_lp):
            lpDict= io.MakeInputDictionary(fname_lp)
        else:
            raise Exception('Error: Long Pulse Input File does not exist')

        # check that the output file exists/make a directory
        if not os.path.exists(oname):
            os.makedirs(oname)


        # SRK - setting an altitude grid
        LPminAlt = self.LP_MinAltitude
        LPmaxAlt = self.LP_MaxAltitude
        tm1=numpy.arange(self.NW_MinAltitude,self.NW_MaxAltitude,self.NW_DeltaAltitude)
        #	tm2=scipy.arange(200.0,400.0,50)
        #	htout=scipy.zeros((tm1.size+tm2.size,2),dtype='float64')
        #	htout[:,0]=scipy.concatenate((tm1,tm2))
        #	htout[:,1]=scipy.concatenate((tm1+tm1[1]-tm1[0],tm2+tm2[1]-tm2[0]))
        htout=numpy.zeros((tm1.size,2),dtype='float64')
        htout[:,0]=tm1
        htout[:,1]=tm1+tm1[1]-tm1[0]
        htoutm=numpy.mean(htout,axis=1)

        outDict = io.MakeOutputDictionary(acDict,htout)

        BMCODES = acDict['BMCODES']
        Ibm = numpy.where((BMCODES[:,2]>=50.0) & (BMCODES[:,2]<=90.0) & (BMCODES[:,1]>=-180.0) & (BMCODES[:,1]<=180.0))[0]
        Nbeams1 = len(Ibm)
        print Nbeams1


        for itime in range(acDict['UnixTime'].shape[0]):
            # getting all the Alternating code

            print 'itime', itime,acDict['UnixTime'].shape[0]

            AllBabs = numpy.ravel(acDict['Babs'][Ibm,:])
            Alldip = numpy.ravel(acDict['dip'][Ibm,:])
            Alldec = numpy.ravel(acDict['dec'][Ibm,:])
            Allk = numpy.reshape(acDict['k1'][Ibm,:,:], (Nbeams1*acDict['NAltitudes'],3))
            AllAltitude = numpy.ravel(acDict['Altitude'][Ibm,:])/1000. # convert to km
            AllVlos = numpy.ravel(acDict['Vlos'][itime,Ibm,:])
            AlldVlos = numpy.ravel(acDict['dVlos'][itime,Ibm,:])
            Allmob = numpy.ravel(acDict['mob'][itime,Ibm,:])
            Allkappa = numpy.ravel(acDict['kappa'][itime,Ibm,:])


            # setting upper limit at 150 km
            Iht = numpy.where(AllAltitude>htout[-1,-1])[0]
            AllVlos[Iht] = numpy.nan

            # grabbing corresponding Long Pulse data
            # repeat the data over the interval
            LPUnixTimeMean = numpy.mean(lpDict['UnixTime'], axis=1)
            I = numpy.where((acDict['UnixTime'][itime,0]<=LPUnixTimeMean) \
                        & (acDict['UnixTime'][itime,1]>=LPUnixTimeMean))[0]
            I = numpy.unique(I)

            # assume the beam pattern is the same between LP and AC
            #repeat this over the

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



            print 'itime', itime
            print 'I', I
            print 'LPabs.shape', LPBabs.shape
            print 'LPAbs', LPBabs.shape
            print 'LPdec', LPdec.shape
            print 'LPdip', LPdip.shape
            print 'LPk1', lpDict['k1'].shape
            print 'LPk1', LPk1.shape, Ibm.shape
            print 'LPk1 2', LPk1.shape, Ibm.shape[0], len(I),lpDict['NAltitudes']
            print 'LPk1', LPk1.shape
            print 'vlos shape', lpDict['Vlos'].shape
            print 'LPVlos', LPVlos.shape
            print 'LPdVlos', LPdVlos.shape
            print 'LPmob', LPmob.shape
            print 'LPkappa', LPkappa.shape
        # io.SaveOutDict(oname,inDict,outDict)

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

            # print 'AllBabs', AllBabs.shape
            # print 'Alldec', Alldec.shape
            # print 'Alldip', Alldip.shape
            # print 'Allk', Allk.shape
            # print 'AllAltitude', AllAltitude.shape
            # print 'AllVlos', AllVlos.shape
            # print 'AlldVlos', AlldVlos.shape
            # print 'Allmob', Allmob.shape
            # print 'Allkappa', Allkappa.shape


            # averaging over the declination and the dip angle
            decAltGrid=numpy.zeros(htout.shape[0],dtype='Float64')
            dipAltGrid=numpy.zeros(htout.shape[0],dtype='Float64')
            BabsAltGrid=numpy.zeros(htout.shape[0],dtype='Float64')
            for bb in range(htout.shape[0]):
                I=numpy.where((AllAltitude>=htout[bb,0])&(AllAltitude<=htout[bb,1]))[0]
                decAltGrid[bb]=numpy.nanmean(Alldec[I])
                dipAltGrid[bb]=numpy.nanmean(Alldip[I])
                BabsAltGrid[bb]=numpy.nanmean(AllBabs[I])


            # invert the winds
            # 09-17-2017 will need to work through this to know what issues are and fix them.
            try:
                (htout,EstimatedWinds,CovEstimatedWind,LosVelEstimated,Iout) = winds.invertwinds(htout,AllAltitude,AllVlos,AlldVlos,Allk,\
                                                     Allmob,Allkappa,decAltGrid,dipAltGrid)
            except:
                print 'failed to run invertWinds'
                break
            # print 'htout', htout
            # print 'EstimatedWinds', EstimatedWinds
            # print 'Covariance Estimated Winds', CovEstimatedWind
            # print 'LOs vel Estimated', LosVelEstimated

            # first 3 elements are the electric field
            # remaining elements are the wind components on the htout grid
            # should be 3+htout.shape[0] elements long
            outDict['Efield'][itime,:] = numpy.transpose(EstimatedWinds[0:3])
            outDict['errEfield'][itime,:] = numpy.sqrt(numpy.diag(EstimatedWinds)[0:3]) # grab the diagonal elements
            tempWind = EstimatedWinds[3:]
            tempErrWind = numpy.sqrt(numpy.diag(CovEstimatedWind)[3:])

            #perp east, north, up
            outDict['WindGmag'][itime,:,0] = tempWind[0::3][numpy.newaxis,:,0]
            outDict['WindGmag'][itime,:,1] = tempWind[1::3][numpy.newaxis,:,0]
            outDict['WindGmag'][itime,:,2] = tempWind[2::3][numpy.newaxis,:,0]
            outDict['errWindGmag'][itime,:,0] = tempErrWind[0::3]
            outDict['errWindGmag'][itime,:,1] = tempErrWind[1::3]
            outDict['errWindGmag'][itime,:,2] = tempErrWind[2::3]

            WindGmag = outDict['WindGmag']
            for ii in range(htout.shape[0]):
                outDict['WindGeo'][itime,ii,:] = numpy.squeeze(winds.gmag2geo(\
                                                numpy.array([[WindGmag[itime,ii,0]],\
                                                [WindGmag[itime,ii,1]],\
                                                [WindGmag[itime,ii,2]]]),\
                                                numpy.deg2rad(decAltGrid[ii]),\
                                                numpy.deg2rad(dipAltGrid[ii])))
            terrgeo = winds.gmag2geo_covar(CovEstimatedWind[3:,3:],\
                                            numpy.deg2rad(decAltGrid),\
                                            numpy.deg2rad(dipAltGrid))
            tempErrWindGeo = numpy.sqrt(numpy.diag(terrgeo))
            outDict['errWindGeo'][itime,:,0] = tempErrWindGeo[0::3]
            outDict['errWindGeo'][itime,:,1] = tempErrWindGeo[1::3]
            outDict['errWindGeo'][itime,:,2] = tempErrWindGeo[2::3]

            # now need to calculate the velocities in the E-region
            # need this to calculate Joule Heating and currents
            VCovarIn = [3000.*3000.,3000.*3000.,30.*30.]
            pIn = [200.0,0.5,2000.0,100.0]
            # PLAT,AllVlos,AlldVlos,Allk,AllPlat,AllPlong,Allht,htmin=150.0*1000,htmax=400.0*1000,\
        	    # covar=[1000.*1000.,1000.*1000.,5.*5.],p=[200.0,0.5,2000.0,100.0])
            (alt_out1,tempVest,tempdVest,tempdVestAll,Foo) = winds.compute_velvec2(htout,\
                                                                    AllVlos,AlldVlos,Allk,\
                                                                    AllAltitude, [],AllAltitude,\
                                                                    htmin=0.0,htmax=1000.0*1000,\
                                                                    covar=VCovarIn,p=pIn)

            # print 'tmpVest.shape', tempVest.shape
            outDict['VestGmag'][itime,:,:] = tempVest
            outDict['errVestGmag'][itime,:,:] = tempdVest

            # note this is different from what Mike N Did
            # you don't which beams were chosen, so you need to loop through the chosen
            # beams and get their index correct

            for ibeam in range(Nbeams1):#Ibm:
                #htoutm = height out mean
                Ve_onNeGrid = interp1d(htoutm,outDict['VestGmag'][itime,:,0],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
                Vn_onNeGrid = interp1d(htoutm,outDict['VestGmag'][itime,:,1],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
                Ue_onNeGrid = interp1d(htoutm,outDict['WindGmag'][itime,:,0],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
                Un_onNeGrid = interp1d(htoutm,outDict['WindGmag'][itime,:,1],bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)
                tempB = interp1d(htoutm,BabsAltGrid,bounds_error=0)(acDict['Altitude'][ibeam,:]/1000.)

                # fill the variables
                # see equation Thayer 1998, JGR, equation 8
                #east = x, north = y
                outDict['CurrentGmag'][itime,ibeam,:,0] = elementary_charge*acDict['Ne'][itime,ibeam,:]*(Ve_onNeGrid + outDict['Efield'][itime,1]/tempB)
                outDict['CurrentGmag'][itime,ibeam,:,1] = elementary_charge*acDict['Ne'][itime,ibeam,:]*(Vn_onNeGrid - outDict['Efield'][itime,0]/tempB)
                # joule heating, electric field only
                outDict['JouleHeatingE'][itime,ibeam,:] = acDict['PedersenConductivity'][itime,ibeam,:]*(outDict['Efield'][itime,:]**2).sum()
                # Joule heating with neutral wind
                # see equation Thayer 1998, JGR, equation 1
                outDict['JouleHeatingTotal'][itime,ibeam,:] = acDict['PedersenConductivity'][itime,ibeam,:]*\
                                                            ((outDict['Efield'][itime,0] - Un_onNeGrid*tempB)**2 + \
                                                            (outDict['Efield'][itime,1] - Ue_onNeGrid*tempB)**2)

                # 09/19/2017 - I am going to hold off on integrating
                #I am going to hold off on doing the errors and the integrated terms
                # only because the integrated terms should really be done along the
                # the field line
                print '\n ibeam', ibeam
                print 'tmpVeast interpolated', Ve_onNeGrid
                print 'Vest East original data', tempVest[:,0]
                # print 'Alt', acDict['Altitude'][ibeam,:]/1000.

            # # print 'errWindGeo', outDict['errWindGeo']
            print 'Ibm', Ibm
        return outDict, acDict

#main(fname_ac,fname_lp,oname)


    # outDict[]

    # Efield[itm,:]=scipy.transpose(test[0:3])
    # errEfield[itm,:]=scipy.sqrt(scipy.diag(terr)[0:3])
    # twind=test[3:]
    # errWind=scipy.sqrt(scipy.diag(terr)[3:])
