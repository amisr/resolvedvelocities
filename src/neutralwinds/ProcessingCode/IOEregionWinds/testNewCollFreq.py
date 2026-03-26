import numpy
import tables

v_elemcharge = 1.602177E-19
v_amu = 1.6605402e-27

def readafile(fname):
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

def compute_ion_neutral_collfreq(densities, Tn, mi, Ti=None):
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
    print sqrtTr.shape, log10Tr.shape
    # pseudocode from ashton
    # inds1 = np.where(Tr > 800)
    # inds2 = np.where(Tr <= 800)
    #
    # nu_in[inds1] += 1.0e-6 * densities[5] * 2.59e-11 * sqrtTr[inds1] * (1.0 - 0.073 * log10Tr[inds1])**2.0
    # nu_in[inds2] += densities[5] * O2const / numpy.sqrt(mi[inds2] * (mi[inds2] + 32.0))
    if mi == 32.0:
        foo1 = numpy.zeros(sqrtTr.shape)
        inds1 = numpy.where(Tr>800)
        # print inds1
        # print  sqrtTr[inds1].shape
        # print log10Tr[inds1].shape
        # print densities[5].shape
        foo1[inds1] = 1.0e-6 * densities[5][inds1] * 2.59e-11 * sqrtTr[inds1] * (1.0 - 0.073 * log10Tr[inds1])**2.0
        print 'foo1', foo1.shape
        print foo1
        nu_in[inds1] += 1.0e-6 * densities[5][inds1] * 2.59e-11 * sqrtTr[inds1] * (1.0 - 0.073 * log10Tr[inds1])**2.0

        inds2 = numpy.where(Tr<=800)
        nu_in[inds2] += densities[5][inds2] * O2const / numpy.sqrt(mi * (mi + 32.0))
        # nu_in[inds2] += densities[5] * O2const / numpy.sqrt(mi[inds2] * (mi[inds2] + 32.0))
    else:
        nu_in += densities[5] * O2const / numpy.sqrt(mi * (mi + 32.0))

    return nu_in



dat1=readafile('/data0/NSF_Winds_ReRunOct2019/2012/06/IPY17/20120615.001/20120615.001_ac_10min-fitcal.h5')

nAr = dat1['/MSIS']['nAr']
nH = dat1['/MSIS']['nH']
nHe = dat1['/MSIS']['nHe']
nN = dat1['/MSIS']['nN']
nN2 = dat1['/MSIS']['nN2']
nNO = dat1['/MSIS']['nNO']
nO = dat1['/MSIS']['nO']
nO2 = dat1['/MSIS']['nO2']
Tn = dat1['/MSIS']['Tn']
mass=dat1['/FittedParams']['IonMass']
fraction = dat1['/FittedParams']['Fits'][:,:,:,0:-1,0]
nuinOrg=dat1['/FittedParams']['Fits'][:,:,:,0:-1,2]
fits = dat1['/FittedParams']['Fits']
Babs1=dat1['/Geomag']['Babs'];
if Babs1[0,0]<1.0e-5:
    Babs1=Babs1*1.0e5 # converting into Telsa?

# mob=self.v_elemcharge/(self.v_amu*nuin); mob[:,:,:,0]=mob[:,:,:,0]/mass[0]; mob[:,:,:,1]=mob[:,:,:,1]/mass[1]
#
#
# kappa=mob*1.0; kappa[:,:,:,0]=kappa[:,:,:,0]*Babs1; kappa[:,:,:,1]=kappa[:,:,:,1]*Babs1
# # mass-weighted params
# nuin1=nuin[:,:,:,0]*fraction[:,:,:,0] + nuin[:,:,:,1]*fraction[:,:,:,1] # mass weighted collision frequency
# mob1=mob[:,:,:,0]*fraction[:,:,:,0]+mob[:,:,:,1]*fraction[:,:,:,1] # mass weighted mobility
# kappa1=kappa[:,:,:,0]*fraction[:,:,:,0]+kappa[:,:,:,1]*fraction[:,:,:,1] # mass weighted kappa


 #H, He, N, O, N2, O2
densities = numpy.array([nH,nHe,nN,nO,nN2,nO2])
nuin = numpy.zeros([mass.shape[0],nO.shape[0],nO.shape[1],nO.shape[2]])
mob = numpy.zeros([mass.shape[0],nO.shape[0],nO.shape[1],nO.shape[2]])
for imass in range(len(mass)):
    print mass[imass]
    tmpnuin = compute_ion_neutral_collfreq(densities,Tn,mass[imass])
    nuin[imass,:] = tmpnuin
    print 'tmpnuin',tmpnuin
    mob[imass,:] = v_elemcharge/(v_amu*tmpnuin)/mass[imass]

kappa = mob*Babs1
mob = numpy.moveaxis(mob,0,-1)
nuin = numpy.moveaxis(nuin,0,-1)
kappa = numpy.moveaxis(kappa,0,-1)


print 'nuin shape', nuin.shape
print 'nuin Org', nuinOrg.shape
print 'kappa shape', kappa.shape
print 'fraction shape', fraction.shape


# now do the summation
nuin1 = numpy.nansum(fraction*nuin, axis=-1)
mob1 = numpy.nansum(fraction*mob, axis=-1)
kappa1 = numpy.nansum(fraction*kappa, axis=-1)
nuinOrg1 = numpy.nansum(fraction*nuinOrg,axis=-1)
# this is how I was originally doing it.
nuin2 = fraction[:,:,:,0]*nuin[:,:,:,0]+fraction[:,:,:,1]*nuin[:,:,:,1]
nuin3 = fraction[:,:,:,0]*nuin[:,:,:,0] + (1.-fraction[:,:,:,0])*nuin[:,:,:,1]
print 'nuin1', nuin1
print 'nuinOrg1', nuinOrg1-nuin1
print 'nuinOrg1 - nuin1', numpy.nanmax(numpy.abs((nuinOrg1-nuin1)/nuinOrg1))

# mob = v_elemcharge/(v_amu*nuin)
