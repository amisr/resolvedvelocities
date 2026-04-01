import numpy
import numpy.linalg

class CalculateEregionWinds:

    '''
    Overarching class used to calculate the neutral winds and with helper parameters
    original translated from process_eregwinds_srk.py
    '''
    def __init__(self):

        return

    def gmag2geo_covar(self,input,dec_all,dip_all,direction=0):

        # Create Jacobian matrix
        Jac=numpy.zeros(input.shape,dtype=input.dtype)

        if direction==0:
            for ii in range(dec_all.size):
                dec=dec_all[ii]
                dip=dip_all[ii]
                Rgmag=numpy.array([
        			[numpy.cos(dec),numpy.sin(dip)*numpy.sin(dec),-numpy.cos(dip)*numpy.sin(dec)],	# east
        			[-numpy.sin(dec),numpy.cos(dec)*numpy.sin(dip),-numpy.cos(dip)*numpy.cos(dec)], # north
        			[0.0,numpy.cos(dip),numpy.sin(dip)]])											# up
                Jac[ii*3:ii*3+3,ii*3:ii*3+3]=Rgmag

        elif direction==1:
            for ii in range(dec_all.size):
                dec=dec_all[ii]
                dip=dip_all[ii]
                Rgeo=numpy.array([
        			[numpy.cos(dec),-numpy.sin(dec),0.0],	# perp east
        			[numpy.sin(dec)*numpy.sin(dip),numpy.cos(dec)*numpy.sin(dip),numpy.cos(dip)], # perp north
        			[-numpy.sin(dec)*numpy.cos(dip),-numpy.cos(dec)*numpy.cos(dip),numpy.sin(dip)]]) # par
                Jac[ii*3:ii*3+3,ii*3:ii*3+3]=Rgeo

        Jac=numpy.matrix(Jac)
        covar=Jac*numpy.matrix(input)*numpy.transpose(Jac)

        return covar

    def gmag2geo(self,input,dec,dip,direction=0):
        # direction:
        # 0 - gmag2geo
        # 1 - geo2gmag

        output=[]
        if direction==0:
            Rgmag=numpy.matrix([
        		[numpy.cos(dec),numpy.sin(dip)*numpy.sin(dec),-numpy.cos(dip)*numpy.sin(dec)],	# east
        		[-numpy.sin(dec),numpy.cos(dec)*numpy.sin(dip),-numpy.cos(dip)*numpy.cos(dec)], # north
        		[0.0,numpy.cos(dip),numpy.sin(dip)]])											# up
            output=Rgmag*numpy.matrix(input)

        elif direction==1:
            Rgeo=numpy.matrix([
        		[numpy.cos(dec),-numpy.sin(dec),0.0],	# perp east
        		[numpy.sin(dec)*numpy.sin(dip),numpy.cos(dec)*numpy.sin(dip),numpy.cos(dip)], # perp north
        		[-numpy.sin(dec)*numpy.cos(dip),-numpy.cos(dec)*numpy.cos(dip),numpy.sin(dip)]]) # par
            output=Rgeo*numpy.matrix(input)

        return output

    def invertwinds(self,htout,ht,Vlos,dVlos,k,mob,kappa,dec,dip,\
        			abserrFilter = 100.,windCovar = [1000.*1000.,1000.*1000,5.*5.],\
        			EfieldCovar = [10.e-1,10.e-1,1.e-5]):

        '''
        This is the main program that is calculating the neutral winds

        Needs to be updated and checked more carefully
        '''

        htmin=htout[0,0]
        htmax=htout[-1,1]

        # print 'vlosin', Vlos

        decall=numpy.zeros(htout.shape[0]+1,dtype='float64')
        dipall=numpy.zeros(htout.shape[0]+1,dtype='float64')
        decall[0]=numpy.nanmean(dec)
        decall[1:]=dec
        dipall[0]=numpy.nanmean(dip)
        dipall[1:]=dip

        decall[numpy.where(numpy.isnan(decall))]=numpy.nanmean(dec)
        dipall[numpy.where(numpy.isnan(dipall))]=numpy.nanmean(dip)

        ht=ht.copy()
        Vlos=Vlos.copy()
        dVlos=dVlos.copy()
        k=k.copy()
        mob=mob.copy()
        kappa=kappa.copy()

        fracerrs=numpy.absolute(dVlos)/(numpy.absolute(Vlos)+abserrFilter)
        abserrs=numpy.absolute(dVlos)

        # get rid of some obvious bad points.
        # this seems to be throwing out the long pulse data?
        # should probably do a consistency check bt >130 km data
        # and data from LP
        I = numpy.where((ht>=htmin) & (ht<=htmax) & (abserrs < abserrFilter) & (numpy.isnan(Vlos) == False))
        # print 'I', I
        ht = ht[I]
        Vlos = Vlos[I]
        dVlos = dVlos[I]
        k = k[I]
        mob = mob[I]
        kappa = kappa[I]
        Iout = I
        print('Vlos after filter', Vlos.shape)

        Nmeas=Vlos.size
        Neqs=htout.shape[0]*3+3

        # initialize output matrices
        xout=numpy.zeros((Neqs),dtype=Vlos.dtype)
        dxout=numpy.zeros((Neqs),dtype=Vlos.dtype)

        # a priori covariance matrix
        windvar = windCovar[0] #1000.0*1000.0
        Zwindvar = windCovar[-1]#5.0*5.0 # if set to 25 made a big difference!
        te = EfieldCovar[0]#[0]10.0e-1
        tep = EfieldCovar[-1]#1e-5
        # create the cov matrix in geographic coordinates
        # 10/03/2017 - need to look into this a bit more - hardcoded stuff?
        # i think this should be redone slightly
        covar=numpy.ones(Neqs,dtype='float64')*windvar
        covar[5::3]=Zwindvar
        covar[0:3]=0.0
        covar[-9:]=Zwindvar # assumes the last 3 altitudes have better covariance - seems hardcoded?
        SigmaV = numpy.array(numpy.diagflat(covar)) # full matrix in geo coords
        # rotate to geomagnetic coordinates
        SigmaV = self.gmag2geo_covar(SigmaV,\
        							numpy.deg2rad(decall),\
        							numpy.deg2rad(dipall),\
        							direction=1)
        # put in electric field guy
        efieldvar=numpy.matrix(numpy.diagflat([te*te,te*te,tep*tep]))
        SigmaV[0:3,0:3]=efieldvar
        # print 'sigmaV.shape 2', SigmaV.shape

        # form C matrix
        C=numpy.array([[1.0/(1.0+kappa*kappa),-kappa/(1.0+kappa*kappa),numpy.zeros(kappa.shape,kappa.dtype)],
        	[kappa/(1.0+kappa*kappa),1.0/(1.0+kappa*kappa),numpy.zeros(kappa.shape,kappa.dtype)],
        	[numpy.zeros(kappa.shape,kappa.dtype),numpy.zeros(kappa.shape,kappa.dtype),numpy.ones(kappa.shape,kappa.dtype)]])

        # build the A matrix
        Amatrix=numpy.matrix(numpy.zeros((Nmeas,Neqs),dtype='float64'))
        for aa in range(Nmeas):
            Amatrix[aa,0:3]=mob[aa]*numpy.matrix(k[aa,:])*numpy.matrix(C[:,:,aa]) # electric field terms
        for aa in range(htout.shape[0]):
            I=numpy.where((ht>=htout[aa,0])&(ht<=htout[aa,1]))[0]
            for bb in range(len(I)):
                Amatrix[I[bb],((aa+1)*3):((aa+1)*3+3)]=numpy.matrix(k[I[bb],:])*numpy.matrix(C[:,:,I[bb]]) # neutral wind terms

        # error covariance matrix
        SigmaE=numpy.matrix(numpy.diagflat(dVlos*dVlos))
        Vlos=numpy.matrix(Vlos[:,numpy.newaxis])
        # print 'Amatrix', Amatrix
        # print 'SigmaV', SigmaV
        # print 'SigmaE', SigmaE
        # print 'Vlos', Vlos
        # print 'Vlos.shape b4 inversion', Vlos.shape
        # print 'Vlos', Vlos

        # do the inversion
        try:
            test = SigmaV*numpy.transpose(Amatrix)*numpy.linalg.inv(Amatrix*SigmaV*numpy.transpose(Amatrix) + SigmaE)*Vlos # estimate
            terr = numpy.linalg.inv(numpy.transpose(Amatrix)*numpy.linalg.inv(SigmaE)*Amatrix + numpy.linalg.inv(SigmaV)) # covariance matrix
            # print 'test', test
        except:
            print('inversion failed in Invert Winds')
            test = numpy.nan*numpy.zeros((Neqs))
            terr = numpy.nan*numpy.zeros((Neqs,Neqs))
        # print 'test', test
        # print 'test shape', test.shape
        # print 'Amatrix shape', Amatrix.shape

        # 09/19/2017 - do not uncomment this line - will get all NaNs on winds and estimated velocities
        VlosEst=numpy.matrix(Amatrix)*numpy.matrix(test) # the forward model
        # print 'VlosEst shape and Altitude.shape', VlosEst.shape, Vlos.shape, ht.shape

        # print 'VlosEst', numpy.ravel(VlosEst), numpy.ravel(Vlos)
        return htout,numpy.array(test),numpy.array(terr),numpy.array(VlosEst),Iout,Vlos,dVlos,ht


    def compute_velvec2(self,PLAT,AllVlos,AlldVlos,Allk,AllPlat,AllPlong,Allht,htmin=150.0*1000,htmax=400.0*1000,\
        covar=[1000.*1000.,1000.*1000.,5.*5.],FracErrorOffset = 200.0, FracErrorThreshold=0.5, AbsoluteErrorThreshold=100.0, \
        p = None ):
        '''
        This is a direct copy of vvels.py and copied from Mike's code directly
    	I cleaned up and made this into numpy functions.

    	Also note, that PLat and altitude can be used synamously.
    	p=[200.0,0.5,2000.0,100.0]
    	'''
    	#
    	# This function takes the line of sight velocities and resolves them into vector velocities as a function of magnetic latitude.
    	# It uses the Bayesian approach and doesn't care about pairs or anything like that.
    	# It does care about the magnetic latitude that you choose for binning, so be wise, eh.
    	#

        if p:
            FracErrorOffset = p[0]
            FracErrorThreshold = p[1]
            AbsoluteErrorThreshold = p[3]

    	# this is the magnetic latitude over which the function bins the data
        plat_in = numpy.copy(PLAT)
        Nplout = plat_in.shape[0]
        plat_out = numpy.zeros((Nplout),dtype='Float64')

        Nparms = Allk.shape[1]
    	# print 'Nparams in vvels'

    	# a priori covariance matrix
        SigmaV = numpy.matrix(numpy.diagflat(covar))

    	# srk comment - Ashton will probably not like this...
        fracerrs = numpy.absolute(AlldVlos)/(numpy.absolute(AllVlos)+FracErrorOffset)
        abserrs = numpy.absolute(AlldVlos)

    	# loop over output latitudes (srk - or altitudes)
        Nmeas = numpy.zeros((Nplout))
        Vest = numpy.zeros((Nplout,Nparms),dtype=AllVlos.dtype)
        dVest = numpy.zeros((Nplout,Nparms),dtype=AllVlos.dtype)
        dVestAll = numpy.zeros((Nplout,Nparms,Nparms),dtype=AllVlos.dtype)
        status = numpy.ones(Nplout)
        for i in range(Nplout):
            plat_out[i] = (plat_in[i,0]+plat_in[i,1])/2.0
    		# print 'plat_out', i, plat_out[i]
            I = numpy.where((AllPlat>plat_in[i,0]) & \
    						(AllPlat<plat_in[i,1]) & \
    						(Allht>htmin) & \
    						(Allht<htmax) & \
    						(numpy.isfinite(AllVlos)) & \
    						(fracerrs <= FracErrorThreshold) & \
    						(abserrs < AbsoluteErrorThreshold))[0]
    		# print 'I'
            Vest[i,:] = numpy.nan
            dVest[i,:] = numpy.nan
            if len(I) != 0:
                try:
                    tvlos = numpy.transpose(numpy.matrix(AllVlos[I]))
                    SigmaE = numpy.matrix(numpy.diagflat(AlldVlos[I]*AlldVlos[I]))
                    A = numpy.matrix(Allk[I,:])
                    tv = SigmaV*numpy.transpose(A)*numpy.linalg.inv(A*SigmaV*numpy.transpose(A) + SigmaE)*tvlos
                    ts = numpy.linalg.inv(numpy.transpose(A)*numpy.linalg.inv(SigmaE)*A + numpy.linalg.inv(SigmaV))
                    Vest[i,:] = numpy.reshape(numpy.array(tv),(1,Nparms))
                    dVest[i,:] = numpy.reshape(numpy.sqrt(numpy.array(numpy.diag(ts))),(1,Nparms))
                    dVestAll[i,:,:]=ts
                    Nmeas[i]=len(I)
                except:
                    status[i] = 0
                    print('failed Vvels Inversion')

        if status.all() == False:
            statusOut = False
        else:
            statusOut = True
    		# print 'len I', len(I), i, Vest[i,:]
    		# else:
    		#	 print 'lenth eq 0'
    		# print 'Vest', Vest

        return plat_out, Vest, dVest, dVestAll, Nmeas,statusOut
