import tables
import copy
import numpy
'''
This class is primarily designed to be the IO into the E-region winds
algorithm.
'''

class IOEregionWinds:
    def __init__(self, nuinScaler = 1.):

        self.v_elemcharge = 1.602177E-19
        self.v_amu = 1.6605402e-27
        self.nuinScaler = nuinScaler
         # set for now will need to change

        self.DictList = {}
        self.DictList['Winds'] = ['WindGmag', 'errWindGmag', 'WindGeo', 'errWindGeo','Altitude']
        self.DictList['VectorVels'] = ['VestGmag','VestGeo','errVestGeo','errVestGmag','Altitude', 'Angle']
        self.DictList['Fregion'] = ['VestGmag_300km','errVestGmag_300km']
        self.DictList['Ne'] = ['MeanSNR', 'MedianSNR', 'MeanNeRaw', 'MedianNeRaw', \
                                'MeanNeFitted', 'MedianNeFitted', 'Altitude', \
                                'NeFittedRaw', 'SNRRaw', 'Kappa', \
                                'VerticalBeamNe', 'errVerticalBeamNe']
        self.DictList['Time'] = ['UnixTime', 'UTDecHrs','LocalDecHrs','MLTDecHrs', 'MeanUnixTimeMe','nuinScaler']
        self.DictList['GeophysicalParameters'] = ['KP','AP','SymH','F107','AE','AL','AU','KPsum','AEmean', 'F107A', 'F107Raw']
        self.DictList['ElectricFields'] = ['Efield','errEfield']
        self.DictList['Ground_Mag'] = ['GroundMag_UnixTime', 'GroundMag_D', 'GroundMag_Z', 'GroundMag_H']
        self.DictList['Status'] = ['Status']
        self.DictList['Forces'] = ['Coriolis', 'Centrifugal', 'Lorentz', 'HallDrag', 'PedersenDrag']
        self.DictList['Information'] = ['Version', 'FileCreationDateTime','ExperimentName','LongPulseFile', \
                                        'AlternatingCodeFile', 'ConfigFile', 'ExperimentDirectory']
        self.DictList['JouleHeating'] = ['VerticalBeamJouleHeatingE', 'VerticalBeamJouleHeatingTotal', \
                                        'MedianJouleHeatingE', 'MedianJouleHeatingTotal', \
                                        'MeanJouleHeatingE', 'MeanJouleHeatingTotal','AltitudeJH', \
                                        'IntegratedMeanJouleHeatingE','IntegratedMeanJouleHeatingTotal',  \
                                        'IntegratedMedianJouleHeatingE','IntegratedMedianJouleHeatingTotal', \
                                        'VerticalBeamJouleHeatingMechanical', \
                                        'IntegratedMeanJouleHeatingMechanical','IntegratedMedianJouleHeatingMechanical', \
                                        'MedianJouleHeatingMechanical', 'MeanJouleHeatingMechanical', \
                                        ]
        self.DictList['JouleHeating_Full'] = ['AltitudeFull', 'JouleHeatingEFull', \
                                                'PedersenConductivtityFull', \
                                                'IntegratedJouleHeatingEFull',\
                                                'PedersenConductanceFull']
        self.DictList['JouleHeating_Thayer'] = ['AltitudeJH', \
                                                'VerticalBeamJouleHeatingTotal_Thayer', \
                                                'VerticalBeamJouleHeatingMechanical_Thayer', \
                                                'VerticalBeamEMTranfer_Thayer', \
                                                'MedianJouleHeatingTotal_Thayer', \
                                                'MeanJouleHeatingTotal_Thayer', \
                                                'MedianJouleHeatingMechanical_Thayer',\
                                                'MeanJouleHeatingMechanical_Thayer', \
                                                'MeanEMTransfer_Thayer', 'MedianEMTransfer_Thayer', \
                                                'IntegratedMedianJouleHeatingTotal_Thayer',
                                                'IntegratedMeanJouleHeatingTotal_Thayer',
                                                'IntegratedMeanJouleHeatingMechanical_Thayer', \
                                                'IntegratedMedianJouleHeatingMechanical_Thayer', \
                                                'IntegratedMeanEMTransfer_Thayer', \
                                                'IntegratedMedianEMTransfer_Thayer'
                                                ]
                                                # 'MedianJouleHeatingTotal_Thayer', 'MedianJouleHeatingMechanical_Thayer', \
                                                # 'MeanJouleHeatingTotal_Thayer', 'MeanJouleHeatingMechanical_Thayer', \
                                                # 'MeanEMTransfer_Thayer', 'MedianEMTransfer_Thayer', \
                                                # 'IntegratedMeanJouleHeatingMechanical_Thayer',\
                                                # 'IntegratedMedianJouleHeatingMechanical_Thayer', \
                                                #
                                                #
                                                #
                                                #
                                                # ]

        self.DictList['Conductivity'] = ['MedianHallConductance', 'MeanHallConductance', \
                                        'MedianPedersenConductance','MeanPedersenConductance', \
                                        'MeanPedersenConductivity', 'MedianPedersenConductivity',
                                        'MeanHallConductivity','MedianHallConductivity', \
                                        'VerticalBeamHallConductivity','VerticalBeamPedersenConductivity', \
                                        'AltitudeJHBeam', 'ScaleHeight', \
                                        'Vertical_nuin', 'Vertical_nuin_Brekke',\
                                        'VerticalTi', 'VerticaldTi', 'VerticalTn']
        self.DictList['ChiSquareTest'] = ['Vlos','dVlos','VlosEst', 'VlosAltGrid', 'dVlosAltGrid']

        return

    def readafile(self,fname):
        '''
        Function written by nicolls to read
        '''
        h5file=tables.open_file(fname)
        output={}
        for group in h5file.walk_groups("/"):
            output[group._v_pathname]={}
            for array in h5file.list_nodes(group, classname = 'Array'):
                output[group._v_pathname][array.name]=array.read()
        h5file.close()

        return output

    def write_outputfile(self,fhandle,dict2do,keys2do=[],groupname='',name=''):
        """
        Writing the outputs
        """
        if groupname == '':
        	group=fhandle.root
        else:
        	if fhandle.__contains__('/'+groupname):
        		group='/'+groupname
        	else:
        		group=fhandle.create_group(fhandle.root, groupname, 'Dataset')

        if len(keys2do)==0:
            try:
                fhandle.removeNode(group,name)
            except:
                ''
            print 'writing output', group, name
            fhandle.create_array(group,name, dict2do, "Dataset")
        else:
            for key in keys2do:
                fhandle.create_array(group, key, dict2do[key], "Dataset")


    def compute_collfreq(self,nO,nN2,nO2,Tn, Ti=None,Te=1000.0, mj=30.0):
        """


        This was taken from model_utils.py as part of the standard ISR Fitter
        SRK made edits to some of the variables to make it more readable


    #   nu_in = 0.0
    #   nu_in = nu_in + d[1]*scipy.sqrt(0.79/16) # O
    #   nu_in = nu_in + d[2]*scipy.sqrt(1.76/28) # N2
    #   nu_in = nu_in + d[3]*scipy.sqrt(1.59/32) # O2
    #   nu_in = nu_in*2.6e-9
    # d[0] - HE NUMBER DENSITY(CM-3)
    #   d[1] - O NUMBER DENSITY(CM-3)
    #   d[2] - N2 NUMBER DENSITY(CM-3)
    #   d[3] - O2 NUMBER DENSITY(CM-3)
        """
        if Ti is not None:
            if Tn.shape[0] == Ti.shape[0]:
                Tr = (Tn+Ti)/2.0

            else:
                raise Exception( 'altitude arrays are not equal for Ti and Tn')
        else:
            print 'in else'
            Ti=Tn
            Tr = (Tn+Ti)/2.0
        nu_in = 0.0
        if mj==16.0: # O
            #O+ + O
            # n(O) = d[1]
            # Schunk and Nagy Table 4.5
            nu_in = nu_in + nO*3.67e-11*numpy.sqrt(Tr)*(1.0-0.064*numpy.log10(Tr))**2.0
        else:
            nu_in = nu_in + 1.0e6*nO*9.14e-15/numpy.sqrt(mj*(mj+16.0)) # not entirely sure where this comes from
        if mj==28.0: # N2
            # n(N2) = d[2]
            #N2+ + N2
            # Schunk and Nagy Table 4.5
            nu_in = nu_in + nN2*5.14e-11*numpy.sqrt(Tr)*(1.0-0.069*numpy.log10(Tr))**2.0
        else:
            nu_in = nu_in + 1.0e6*nN2*1.80e-14/numpy.sqrt(mj*(mj+28.0))

        if mj==32.0:# and Tn>800.0: # O2
            #n(O2) = d[3]
            # O2+ + O2
            # Schunk and Nagy Table 4.5
            nu_in = nu_in + nO2*2.59e-11*numpy.sqrt(Tr)*(1.0-0.073*numpy.log10(Tr))**2.0
        else:
            nu_in = nu_in + 1.0e6*nO2*1.83e-14/numpy.sqrt(mj*(mj+32.0))

        nu_en = 0.0
        nu_en = nu_en +nO*8.2e-10*numpy.sqrt(Te) # O
        nu_en = nu_en + nN2*2.33e-11*(1-1.2e-4*Te)*Te # N2
        nu_en = nu_en + nO2*1.8e-10*(1+3.6e-2*numpy.sqrt(Te))*numpy.sqrt(Te)

        return nu_in, nu_en


    """
    Copied from https://github.com/amisr/flipchem/blob/master/flipchem/msis.py
    on 7 April 2020
    """
    def compute_ion_neutral_collfreq(self,densities, Tn, mi, Ti=None):
        """This code calculates the elastic and resonant ion-neutral collision
        frequencies following Chapter 4 of [1]_.
        Parameters
        ==========
        densities : array_like
            An array of neutral densities in this order: H, He, N, O, N2, O2 with
            units of number per cubic meter.
        Tn : float
            The mean neutral temperature in Kelvin
        mi : integer
            The ion mass in amu
        Ti : float, optional
            The ion temperature in Kelvin
        Returns
        =======
        nu_in : float
            The the total ion-neutral collision frequency summed over the
            collisions between the input ion and H, He, N, O, N2, O2
        References
        ==========
        .. [1] Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
               Physics, and Chemistry (Cambridge Atmospheric and Space Science
               Series). Cambridge: Cambridge University Press. 99-054707
               ISBN: 0 521 60770 1
        """
        # set the ion temperature if None
        if Ti is None:
            Ti = Tn

        # Resonant collision frequencies (mi==##) come from
        # table 4.5 page 99 of Schunk and Nagy Ionospheres text book 2000
        # Non-resonant collisions come from equation 4.88 on page 83 and
        # table 4.1. There's some unit conversion needed, but basically
        # the equations come from:
        #    n_in = 2.21*pi*n_n*sqrt(gamma_n*e**2*m_n)/sqrt(m_i*(m_i+m_n))
        # which reduces to (after unit conversion to SI):
        #    n_in = 2.58790619679528e-3*n_n*sqrt(gamma_n*m_n)/sqrt(m_i*(m_i+m_n))
        # where gamma_n is in units of 1e-24 cm^3 and the masses are in amu
        # then, taking the values for gamma in table 4.1 and the amu for the
        # neutrals, we can get an equation that is:
        #    n_in = const * n_n / sqrt(m_i*(m_i+m_n))
        # for each ion. Here's the constants for H+, He+, N+, O+, N2+, and O2+:
        Hconst  = 2.118436361550408e-15
        Heconst = 2.372016271579909e-15 # assuming 100% He-4
        Nconst  = 1.0293931179488746e-14
        Oconst  = 9.0841307137989e-15
        N2const = 1.8168261427597804e-14
        O2const = 1.8518806819759527e-14

        # define some other constants
        Tr = (Tn + Ti) / 2.0
        sqrtTr = numpy.sqrt(Tr)
        log10Tr = numpy.log10(Tr)

        # Now calculate the total ion-neutral collision frequency
        nu_in = 0.0

        # H, no resonant because we don't include H+ in ISR fitting
        nu_in += densities[0] * Hconst / numpy.sqrt(mi * (mi + 16.0))
        # He, no resonant because we don't include He+ in ISR fitting
        nu_in += densities[1] * Heconst / numpy.sqrt(mi * (mi + 16.0))
        # N
        if mi == 14.0:
            nu_in += 1.0e-6 * densities[2] * 3.83e-11 * sqrtTr * (1.0 - 0.063 * log10Tr)**2.0
        else:
            nu_in += densities[2]* Nconst / numpy.sqrt(mi * (mi + 16.0))
        # O
        if mi == 16.0:
            nu_in += 1.0e-6 * densities[3] * 3.67e-11 * sqrtTr * (1.0 - 0.064 * log10Tr)**2.0
        else:
            nu_in += densities[3] * Oconst / numpy.sqrt(mi * (mi + 16.0))
        # N2
        if mi == 28.0:
            nu_in += 1.0e-6 * densities[4] * 5.14e-11 * sqrtTr * (1.0 - 0.069 * log10Tr)**2.0
        else:
            nu_in += densities[4] * N2const / numpy.sqrt(mi * (mi + 28.0))
        # 02
        # SRK note, I changed Tr to Tr.all()
        if mi == 32.0:
             # and Tr.all() > 800.0:
            #srk comment out
            inds1 = numpy.where(Tr>800)
            nu_in[inds1] += 1.0e-6 * densities[5][inds1] * 2.59e-11 * sqrtTr[inds1] * (1.0 - 0.073 * log10Tr[inds1])**2.0
            #nu_in += 1.0e-6 * densities[5] * 2.59e-11 * sqrtTr * (1.0 - 0.073 * log10Tr)**2.0

            inds2 = numpy.where(Tr<=800)
            nu_in[inds2] += densities[5][inds2] * O2const / numpy.sqrt(mi * (mi + 32.0))
        else:
            nu_in += densities[5] * O2const / numpy.sqrt(mi * (mi + 32.0))

        return nu_in




    def compute_electron_neutral_collfreq(self,densities, Te):
        """
        This code calculates the elastic electron-neutral collision frequencies
        following Chapter 4 of [2]_.
        Parameters
        ==========
        densities : array_like
            An array of neutral densities in this order: H, He, N, O, N2, O2 with
            units of number per cubic meter.
            See Notes for comment about Nitrogen.
        Te : float
            The electron temperature in Kelvin
        Returns
        =======
        nu_en : float
            The the total electron-neutral collision frequency summed over the
            collisions between electrons and H, He, O, N2, O2.
        Notes
        =====
        The current code DOES NOT include an electron-Nitrogen collision frequency.
        References
        ==========
        .. [2] Schunk, R., & Nagy, A. (2000). Ionospheres: Physics, Plasma
               Physics, and Chemistry (Cambridge Atmospheric and Space Science
               Series). Cambridge: Cambridge University Press. 99-054707
               ISBN: 0 521 60770 1
        """

        # Elastic electron-neutral collision frequencies only
        # Table 4.6 page 99 of Schunk and Nagy Ionospheres text book 2000
        sqrtTe = numpy.sqrt(Te)
        nu_en = 0.0
        nu_en += 1e-6 * densities[0] * 4.5e-9 * (1.0 - 1.35e-4 * Te) * sqrtTe # H
        nu_en += 1e-6 * densities[1] * 4.6e-10 * sqrtTe # He
        # no N in Schunk and Nagy!
        nu_en += 1e-6 * densities[3] * 8.9e-11 * (1.0 + 5.7e-4 * Te) * sqrtTe # O
        nu_en += 1e-6 * densities[4] * 2.33e-11 * (1.0 - 1.2e-4 * Te) * Te # N2
        nu_en += 1e-6 * densities[5] * 1.82e-10 * (1.0 + 3.6e-2 * sqrtTe) * sqrtTe # O2

        return nu_en



    def ProcessConductivity(self,fname_ac, opt='Original'):
        """
        Added on 11/14/2018 from FastConductivity.py
        passing in the h5 pointer
        """
        # bring this up to speed with python2.7
        # read in the data
        with tables.open_file(fname_ac) as h5:


            ne1=h5.root.FittedParams.Ne.read()#dat1['/FittedParams']['Ne']
            dne1=h5.root.FittedParams.dNe.read()# dne1=dat1['/FittedParams']['dNe']
            time1=h5.root.Time.UnixTime.read()#dat1['/Time']['UnixTime']
            dtime1=h5.root.Time.dtime.read()#dat1['/Time']['dtime']
        # (Nrecs1,Nbeams1,Nhts1)=vlos1.shape
            Altitude=h5.root.FittedParams.Altitude.read()#dat1['/FittedParams']['Altitude']

            Babs1=h5.root.Geomag.Babs.read()#dat1['/Geomag']['Babs'];
            if Babs1[0,0]<1.0e-5:
                Babs1=Babs1*1.0e5 # convert to nanoTelsa
        # '/FittedParams/IonMass' : [('TITLE','Ion Mass'),('Description', 'Mass of ions used in fitting'),('Unit','amu')],\
        # '/FittedParams/Fits' : [('TITLE','Fits'),('Description','Fitted parameters'),('Size','Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, collision frequency, LOS speed)'),('Unit','N/A, Kelvin, s^{-1}, m/s')],\
            mass=h5.root.FittedParams.IonMass.read()#dat1['/FittedParams']['IonMass']
            Fits=h5.root.FittedParams.Fits.read()#dat1['/FittedParams']['Fits'][:,:,:,0:2,0] # O+ fraction
            BeamCodes = h5.root.BeamCodes.read()
            # convert to cm^-3
            nO = h5.root.MSIS.nO.read()/1e6
            nO2 = h5.root.MSIS.nO2.read()/1e6
            nN2 = h5.root.MSIS.nN2.read()/1e6
            Tn = h5.root.MSIS.Tn.read()
            IonMass = h5.root.FittedParams.IonMass.read()
        #nuin=dat1['/FittedParams']['Fits'][:,:,:,0:2,2]

        fraction = Fits[:,:,:,0:-1,0] # fraction including electrons as index = -1
        nuin = Fits[:,:,:,0:-1,2] # collision frequency, including electrons as index = -1
        nuen = Fits[:,:,:,-1,2]
        Ti = Fits[:,:,:,0:-1,1]
        Te = Fits[:,:,:,-1,1]
        PedCond = numpy.zeros(ne1.shape)
        HallCond = numpy.zeros(ne1.shape)

        v_elemcharge = self.v_elemcharge
        v_amu = self.v_amu


        # print nO.shape, nN2.shape, nO2.shape, Tn.shape, Ti[:,:,:,0].shape
        # print 'fraction', fraction
        nuinISR = numpy.zeros(Ti.shape)*numpy.nan


        for i in range(mass.shape[0]):
            print mass[i]
            tmpnuin,tmpnuen = self.compute_collfreq(nO,nN2,nO2,Tn,Ti=Ti[:,:,:,i],mj=mass[i])
            # print 'tmp nu in',tmpnuin.shape
            nuinISR[:,:,:,i] = tmpnuin


        if opt == 'Original':
            mob=v_elemcharge/(v_amu*nuin); # original form in Keely textbook 2009 p. 42
        if opt == 'ISR':
            mob = v_elemcharge/(v_amu*nuinISR)
        # nuinScaler = 1.0 # scalar to increase value of ion-neutral collision frequency
        Babs = numpy.tile(Babs1,ne1.shape[0]).reshape(ne1.shape) # generate a time x beam x altitude array of magnetic field


        # loop over the ion type.
        # assume that since nuIn is weighted by fractional amount that
        # performs any weighting in the conductivity equation
        for i in range(mass.shape[0]):
            mob[:,:,:,i] = mob[:,:,:,i]/mass[i] # equivalent to the form on p 42 of Kelley 2009 textbook
            kappa = mob*1.0
            kappa[:,:,:,i] = mob[:,:,:,i]*Babs1 # equivalent to the form on p 42 of Kelley 2009 textbook

            # equation 2.40a in Kelley's textbook
            # note that the weighting is on the outside equivalent to Evans 1977, JGR
            # weighting is fractional concentration, so from 0-1.
            sp1 = fraction[:,:,:,i]*ne1*v_elemcharge*v_elemcharge/(v_amu*mass[i]*nuin[:,:,:,i]*(1.0+(kappa[:,:,:,i])**2.0))

            # equation 10 from Bostrom 1964 in my notebook
            # should be good > 85 km
            # note again the weighting on the outside
            sh1 = fraction[:,:,:,i]*ne1*v_elemcharge/(Babs*(1.0+(kappa[:,:,:,i])**2.0)) # fraction in effect weights things
            # print sp1
            PedCond = PedCond + sp1
            HallCond = HallCond + sh1

        return PedCond, HallCond





    def MakeInputDictionary(self,fname):
        '''
        This routine is designed to read in the data and reform it into a dictionary

        '''

        inDict = {}
        dat1=self.readafile(fname)

    	# alternating code
    	# '/FittedParams/Ne' : [('TITLE','Electron Density'),
        #('Description', 'Fitted electron density'),
        #('Size','Nrecords x Nbeams x Nranges'),('Unit','m^{-3}')],\
         # '/FittedParams/Fits' : [('TITLE','Fits'),('Description','Fitted parameters'),
         # ('Size','Nrecords x Nbeams x Nranges x Nions+1 x 4 (fraction, temperature, collision frequency, LOS speed)'),('Unit','N/A, Kelvin, s^{-1}, m/s')],\
        # fitted Parameters information
        inDict['Ne'] = dat1['/FittedParams']['Ne']
        ne1 = dat1['/FittedParams']['Ne']
    	inDict['dNe'] = dat1['/FittedParams']['dNe']
    	inDict['Vlos']= dat1['/FittedParams']['Fits'][:,:,:,0,3]
        vlos1 = dat1['/FittedParams']['Fits'][:,:,:,0,3]
        inDict['Ti'] = dat1['/FittedParams']['Fits'][:,:,:,0,1]
        inDict['dTi'] = dat1['/FittedParams']['Errors'][:,:,:,0,1]
    	inDict['dVlos']= dat1['/FittedParams']['Errors'][:,:,:,0,3]
        inDict['Altitude']=dat1['/FittedParams']['Altitude']

        # print dat1['/FittedParams']['Fits'][:,:,:,0,3],


        print dat1['/FittedParams/FitInfo'].keys()
        # added fit information 10-15-2019
        inDict['chi2'] = dat1['/FittedParams/FitInfo']['chi2']
        inDict['fitcode'] = dat1['/FittedParams/FitInfo']['fitcode']

        # time information
        inDict['UnixTime']= dat1['/Time']['UnixTime']
    	inDict['dtime']= dat1['/Time']['dtime']
    	inDict['MLT'] = dat1['/Time']['MagneticLocalTimeSite']

        (Nrecs1,Nbeams1,Nhts1)=vlos1.shape
        inDict['Nrecs'] = Nrecs1
        inDict['Nbeams'] = Nbeams1
        inDict['NAltitudes'] = Nhts1

        inDict['NeRaw'] = dat1['/NeFromPower']['Ne_Mod']
        inDict['AltitudeRaw'] = dat1['/NeFromPower']['Altitude']
        inDict['SNRRaw'] = dat1['/NeFromPower']['SNR']


    	kpn1=dat1['/Geomag']['kpn']; kpe1=dat1['/Geomag']['kpe']; kpar1=dat1['/Geomag']['kpar']
    	k1=numpy.zeros((Nbeams1,Nhts1,3),dtype=kpn1.dtype)
    	k1[:,:,0]=kpe1; k1[:,:,1]=kpn1; k1[:,:,2]=kpar1
        inDict['k1'] = k1
        inDict['kpn'] = kpn1
        inDict['kpe'] = kpe1
        inDict['kpap'] = kpar1

        Babs1=dat1['/Geomag']['Babs'];
    	if Babs1[0,0]<1.0e-5:
    		Babs1=Babs1*1.0e5 # converting into Telsa?

        inDict['Babs'] = Babs1
        inDict['dec'] = dat1['/Geomag']['Declination'];
        inDict['dip'] = dat1['/Geomag']['Dip']

        inDict['nMass'] = dat1['/MSIS']['nMass']
        nMass = dat1['/MSIS']['nMass'] # kg/m^{3}
        nAr = dat1['/MSIS']['nAr']
        nH = dat1['/MSIS']['nH']
        nHe = dat1['/MSIS']['nHe']
        nN = dat1['/MSIS']['nN']
        nN2 = dat1['/MSIS']['nN2']
        nNO = dat1['/MSIS']['nNO']
        nO = dat1['/MSIS']['nO']
        nO2 = dat1['/MSIS']['nO2']
        Tn = dat1['/MSIS']['Tn']
        inDict['Tn'] = Tn
        nTot = nAr+nH+nHe+nN+nN2+nNO+nO+nO2
        totMass = nMass/nTot
        ScaleHeight = 1.38e-23*Tn/totMass/9.8
        inDict['ScaleHeight'] = ScaleHeight

        inDict['F107A'] = dat1['/MSIS']['f107a']
        inDict['F107Raw'] = dat1['/MSIS']['f107']

        mass=dat1['/FittedParams']['IonMass']
        fraction = dat1['/FittedParams']['Fits'][:,:,:,0:-1,0] # O+ fraction

        """
        original way of doing this
        commented out on 10 April 2020

    	#fraction=dat1['/FittedParams']['Fits'][:,:,:,0:3,0] # O+ fraction
    	nuin=dat1['/FittedParams']['Fits'][:,:,:,0:-1,2]

    	mob=self.v_elemcharge/(self.v_amu*nuin);
        mob[:,:,:,0]=mob[:,:,:,0]/mass[0];
        mob[:,:,:,1]=mob[:,:,:,1]/mass[1]
        mob[:,:,:,2]=mob[:,:,:,2]/mass[2]


    	kappa=mob*1.0; kappa[:,:,:,0]=kappa[:,:,:,0]*Babs1; kappa[:,:,:,1]=kappa[:,:,:,1]*Babs1
        kappa[:,:,:,2]=kappa[:,:,:,2]*Babs1

        # mass-weighted params
    	nuin1=nuin[:,:,:,0]*fraction[:,:,:,0] + nuin[:,:,:,1]*fraction[:,:,:,1]+nuin[:,:,:,2]*fraction[:,:,:,2] # mass weighted collision frequency
    	mob1=mob[:,:,:,0]*fraction[:,:,:,0]+mob[:,:,:,1]*fraction[:,:,:,1]+mob[:,:,:,2]*fraction[:,:,:,2] # mass weighted mobility
    	kappa1=kappa[:,:,:,0]*fraction[:,:,:,0]+kappa[:,:,:,1]*fraction[:,:,:,1]++kappa[:,:,:,2]*fraction[:,:,:,2] # mass weighted kappa

        inDict['kappa'] = kappa1
        inDict['mob'] = mob1
        inDict['nuin'] = nuin1


        """
        densities = numpy.array([nH,nHe,nN,nO,nN2,nO2])
        nuin = numpy.zeros([mass.shape[0],nO.shape[0],nO.shape[1],nO.shape[2]])
        mob = numpy.zeros([mass.shape[0],nO.shape[0],nO.shape[1],nO.shape[2]])
        for imass in range(len(mass)):
            print mass[imass]
            tmpnuin = self.compute_ion_neutral_collfreq(densities,Tn,mass[imass])
            tmpnuin = tmpnuin*self.nuinScaler
            nuin[imass,:] = tmpnuin
            # print 'tmpnuin',tmpnuin
            mob[imass,:] = self.v_elemcharge/(self.v_amu*tmpnuin)/mass[imass]

        kappa = mob*Babs1
        mob = numpy.moveaxis(mob,0,-1)
        nuin = numpy.moveaxis(nuin,0,-1)
        kappa = numpy.moveaxis(kappa,0,-1)



        # print 'nuin shape', nuin.shape
        # print 'nuin Org', nuinOrg.shape
        # print 'kappa shape', kappa.shape
        # print 'fraction shape', fraction.shape

        """
        Commented out on 15 April 2020
        # index 0 is O+, index 1 is O2+ and index 2 is NO+
        # No ion neutral collision frequencies for NO+, so what could do instead
        # is to take O+ fraction and take 1-0+fraction to
        # did this as a check from previous work.
        nuin1=nuin[:,:,:,0]*fraction[:,:,:,0] + nuin[:,:,:,1]*fraction[:,:,:,1] # mass weighted collision frequency
    	mob1=mob[:,:,:,0]*fraction[:,:,:,0]+mob[:,:,:,1]*fraction[:,:,:,1] # mass weighted mobility
    	kappa1=kappa[:,:,:,0]*fraction[:,:,:,0]+kappa[:,:,:,1]*fraction[:,:,:,1] # mass weighted kappa


        inDict['kappa'] = kappa1
        inDict['mob'] = mob1
        inDict['nuin'] = nuin1
        """

        # now do the summation
        nuin1 = numpy.nansum(fraction*nuin, axis=-1)
        mob1 = numpy.nansum(fraction*mob, axis=-1)
        kappa1 = numpy.nansum(fraction*kappa, axis=-1)

        inDict['kappa'] = kappa1
        inDict['mob'] = mob1
        inDict['nuin'] = nuin1

        # brekke 1994 equation 5 in units of 1/m^{3}
        nuinBrekke = 4.34e-16*nN2+4.28e-16*nO2+2.44e-16*nO # brekke formula,
        inDict['nuinBrekke'] = nuinBrekke
        inDict['nuinScaler'] = numpy.ones(Nrecs1)*self.nuinScaler




        # commented out on 4/7/2020
        # nuin2 = 4.34e-16*nN2+4.28e-16*nO2+2.44e-16*nO # brekke formula
        # nuin2 = nuin2*3. # need to comment this out in the future
        # meanMass = numpy.mean(mass)
        # mob2 = self.v_elemcharge/(self.v_amu*nuin2*numpy.mean(mass))
        # kappa2 = mob2*Babs1
        #
        # inDict['kappa'] = kappa2
        # inDict['mob'] = mob2
        # inDict['nuin'] = nuin2

        # note not using kappa1
        # calculating the Hall and Pedersen conductivity
        #11 14 2018 original formula commented out
        # sp1=ne1*self.v_elemcharge*self.v_elemcharge/(self.v_amu*mass[0]*self.nuinScaler*nuin[:,:,:,0]*(1.0+(kappa[:,:,:,0]/self.nuinScaler)**2.0))*fraction[:,:,:,0]
    	# +ne1*self.v_elemcharge*self.v_elemcharge/(self.v_amu*mass[1]*self.nuinScaler*nuin[:,:,:,1]*(1.0+(kappa[:,:,:,1]/self.nuinScaler)**2.0))*fraction[:,:,:,1] # Pedersen conductance

        # formula in from fast conductivity.py, which I should probably use that repo altogether.
        try:
            tmpPedCond,tmpHallCond = self.ProcessConductivity(fname, opt='Original')
            inDict['PedersenConductivity'] = tmpPedCond
            inDict['HallConductivity'] = tmpHallCond
        except:
            inDict['PedersenConductivity'] =  dat1['/FittedParams']['Ne']*numpy.nan
            inDict['HallConductivity'] =  dat1['/FittedParams']['Ne']*numpy.nan
            pass

        # will set this with config files
    	BMCODES=dat1['/']['BeamCodes']; Nbeams0=BMCODES.shape[0]
    	Ibm=numpy.where((BMCODES[:,2]>=50.0) & (BMCODES[:,2]<=90.0) & (BMCODES[:,1]>=-180.0) & (BMCODES[:,1]<=180.0))[0]
    	Nbeams1=len(Ibm)
    	# print Nbeams1
    	# print BMCODES[Ibm,:]

        inDict['BMCODES'] = BMCODES
        inDict['Ibm'] = Ibm
        inDict['Nbeams'] = Nbeams1

        # added on 13 March 2018 - need longitude
        inDict['SiteLongitude'] = dat1['/Site']['Longitude']
        inDict['SiteLatitude'] = dat1['/Site']['Latitude']


        return inDict


    def MakeOutputDictionary(self,inDict,AltitudeArr):
        '''
        This makes the output dictionary
        Basically need the input dictionary since it gives the shapes
        of items that will be in the output dictionary
        '''
        outDict = {}

        TimeOnly = numpy.zeros((inDict['UnixTime'].shape[0]), dtype='float64')*numpy.nan
        TimeOnlyList = ['IntegratedVerticalBeamJouleHeatingE', \
                        'IntegratedVerticalBeamJouleHeatingTotal', \
                        'IntegratedMedianJouleHeatingE', \
                        'IntegratedMedianJouleHeatingTotal', \
                        'IntegratedMeanJouleHeatingE', \
                        'IntegratedMeanJouleHeatingTotal', \
                        'IntegratedMeanJouleHeatingMechanical', \
                        'IntegratedMedianJouleHeatingMechanical', \
                        'IntegratedVerticalBeamJouleHeatingMechanical', \
                        'MedianHallConductance', 'MeanHallConductance', \
                        'MedianPedersenConductance','MeanPedersenConductance',\
                        'nuinScaler', 'F107A', 'F107Raw', \
                        'IntegratedMedianJouleHeatingTotal_Thayer',
                        'IntegratedMeanJouleHeatingTotal_Thayer',
                        'IntegratedMeanJouleHeatingMechanical_Thayer', \
                        'IntegratedMedianJouleHeatingMechanical_Thayer', \
                        'IntegratedMeanEMTransfer_Thayer', \
                        'IntegratedMedianEMTransfer_Thayer', \
                        'IntegratedJouleHeatingEFull',\
                        'PedersenConductanceFull']
        for ii in TimeOnlyList:
            outDict[ii] = copy.copy(TimeOnly)


        Timeby3 = numpy.zeros((inDict['UnixTime'].shape[0],3), dtype='float64')*numpy.nan
        Timeby3List = ['Efield', 'errEfield', 'VestGmag_300km','errVestGmag_300km']
        for ii in Timeby3List:
            outDict[ii] = copy.copy(Timeby3)

        # time x altitude x 3 - NW grid
        TimebyAltitudeby3 = numpy.zeros((inDict['UnixTime'].shape[0],AltitudeArr.shape[0],3),\
                            dtype='float64')*numpy.nan
        TimebyAltitudeby3List = ['WindGmag', 'errWindGmag', 'WindGeo', 'errWindGeo', \
                                'VestGmag', 'errVestGmag','VestGeo','errVestGeo', \
                                'Coriolis', 'Centrifugal', 'Lorentz']
        for ii in TimebyAltitudeby3List:
            outDict[ii] = copy.copy(TimebyAltitudeby3)


        NeShape2 = numpy.zeros((inDict['Ne'].shape[0],inDict['Ne'].shape[1],inDict['Ne'].shape[2],2), \
                    dtype='float64')*numpy.nan
        NeShape2List = ['CurrentGmag', 'errCurrentGmag']

        for ii in NeShape2List:
            outDict[ii] = copy.copy(NeShape2)


        # time x altitude - NeutralWind grid
        # wanted to put Joule altitude on same grid
        TimebyAltitude = numpy.zeros((inDict['UnixTime'].shape[0],AltitudeArr.shape[0]),\
                            dtype='float64')*numpy.nan
        TimebyAltitudeList = ['MedianJouleHeatingE', 'MedianJouleHeatingTotal', \
                                'MeanJouleHeatingE', 'MeanJouleHeatingTotal', \
                                'MeanSNR', 'MedianSNR', 'MeanNeRaw', 'MedianNeRaw', \
                                'MeanNeFitted', 'MedianNeFitted', 'VerticalBeamJouleHeatingE', \
                                'VerticalBeamJouleHeatingTotal','ScaleHeight', 'Kappa', \
                                'MedianJouleHeatingMechanical', 'MeanJouleHeatingMechanical', \
                                'VerticalBeamJouleHeatingMechanical', 'MeanPedersenConductivity', \
                                'MedianPedersenConductivity', 'MeanHallConductivity','MedianHallConductivity', \
                                'VerticalBeamHallConductivity','VerticalBeamPedersenConductivity', \
                                'MeanIonNeutralCollisionFrequency', 'StdIonNeutralCollisionFrequency', \
                                'Vertical_nuin', 'Vertical_nuin_Brekke',\
                                'VerticalTi', 'VerticaldTi', 'VerticalTn', 'VerticalBeamNe', 'errVerticalBeamNe', \
                                'PedersenDrag', 'HallDrag', 'Angle', \
                                'VerticalBeamJouleHeatingTotal_Thayer','VerticalBeamJouleHeatingMechanical_Thayer', \
                                'VerticalBeamEMTranfer_Thayer' , 'MedianJouleHeatingTotal_Thayer', \
                                'MeanJouleHeatingTotal_Thayer', \
                                'MedianJouleHeatingMechanical_Thayer',\
                                'MeanJouleHeatingMechanical_Thayer', \
                                'MeanEMTransfer_Thayer', 'MedianEMTransfer_Thayer']
                                #                 'VerticalBeamJouleHeatingTotal', 'VerticalBeamerrJouleHeatingTotal']]
        for ii in TimebyAltitudeList:
            outDict[ii] = copy.copy(TimebyAltitude)

        # time x altitude x 1000 - Raw Ne/Joule Heating Grid
        # added 09/08/2018

        TimebyAltitude = numpy.zeros((inDict['UnixTime'].shape[0],AltitudeArr.shape[0],300),\
                            dtype='float64')*numpy.nan
        TimebyAltitudeList = ['NeFittedRaw', 'SNRRaw','VlosAltGrid', 'dVlosAltGrid']
        for ii in TimebyAltitudeList:
            outDict[ii] = copy.copy(TimebyAltitude)

        # time by 250
        Timeby250 = numpy.zeros((inDict['UnixTime'].shape[0],500),\
                            dtype='float64')*numpy.nan
        Timeby250List = ['Vlos','dVlos','VlosEst']
        for ii in Timeby250List:
            outDict[ii] = copy.copy(Timeby250)



        # time x beam x altitude x 2 - Ne grid
        NeShape2 = numpy.zeros((inDict['Ne'].shape[0],inDict['Ne'].shape[1],inDict['Ne'].shape[2],2), \
                    dtype='float64')*numpy.nan
        NeShape2List = ['CurrentGmag', 'errCurrentGmag']
        for ii in NeShape2List:
            outDict[ii] = copy.copy(NeShape2)


        # time x altitude x beam - Ne grid

        NeShape = numpy.zeros((inDict['Ne'].shape[0],inDict['Ne'].shape[1],inDict['Ne'].shape[2]), \
                  dtype='float64')*numpy.nan
        NeShapeList = ['JouleHeatingE', 'errJouleHeatingE', \
                        'JouleHeatingTotal', 'errJouleHeatingTotal', \
                        'PedersenConductivity', 'HallConductivity', \
                        'JouleHeatingMechanical', 'JouleHeatingTotalThayer', \
                        'JouleHeatingMechanicalThayer', 'EMTransferRateThayer']

        for ii in NeShapeList:
            outDict[ii] = copy.copy(NeShape)


        # time x altitude - Ne grid
        NeShape = numpy.zeros((inDict['Ne'].shape[0],inDict['Ne'].shape[2]), \
                  dtype='float64')*numpy.nan
        # NeShapeList = ['VerticalBeamJouleHeatingE', 'VerticalBeamerrJouleHeatingE', \
        #                 'VerticalBeamJouleHeatingTotal', 'VerticalBeamerrJouleHeatingTotal']
        #                 # can add weighted mean once I have proper uncertainty
        # for ii in NeShapeList:
        #     outDict[ii] = copy.copy(NeShape)

        # time x beam
        NeShape = numpy.zeros((inDict['Ne'].shape[0],inDict['Ne'].shape[1]), \
                  dtype='float64')*numpy.nan
        NeShapeList = ['IntegratedJouleHeatingE', 'IntegratedJouleHeatingTotal',\
                        'IntegratedJouleHeatingMechanical']

        for ii in NeShapeList:
            outDict[ii] = copy.copy(NeShape)

        # added on 5/3/2022
        # time x interpolated beam
        tmpxx = numpy.arange(90.,255.,5.)
        NeShape = numpy.zeros((inDict['Ne'].shape[0],tmpxx.shape[0]), \
                  dtype='float64')*numpy.nan
        NeShapeList = ['JouleHeatingEFull', 'PedersenConductivtityFull']

        for ii in NeShapeList:
            outDict[ii] = copy.copy(NeShape)



        # empty stuff
        NeShapeList = ['GroundMag_UnixTime', 'GroundMag_H', \
                        'GroundMag_D', 'GroundMag_Z', 'Nu', 'F107', \
                        'AE', 'KP', 'AP','KPsum','APmean',\
                        'AL','AU','SymH', 'AltitudeJHBeam', 'AltitudeJH', 'AltitudeFull']

        for ii in NeShapeList:
            outDict[ii] = {}




        return outDict

    def SaveOutDict(self,oname,inDict,outDict):

        """
        Save output dictionary
        """
        outh5file=tables.open_file(oname, mode = "w", title = "Fit File")

        # this contains all of the info you might need to do this

        print 'inDict keys', inDict.keys()
        print 'outDict keys', outDict.keys()

        for ikeys in self.DictList.keys():

            for ilst in self.DictList[ikeys]:

                # 03-14-2018 making a decision that all parameters must be in the outdict...
                # this is potentially confusing what is happening and could result
                # in weird things happening.

                # if ilst in inDict.keys():
                #     print 'ikeys,ilst,shape',ikeys,ilst
                    #self.write_outputfile(outh5file,inDict[ilst],groupname=ikeys,name=ilst)
                if ilst in outDict.keys():
                    print 'ikeys,ilst,shape',ikeys,ilst
                    try:
                        self.write_outputfile(outh5file,outDict[ilst],groupname=ikeys,name=ilst)
                    except:
                        pass

        outh5file.close()
        print 'output location final', oname
        return
