#!/usr/bin/env python

"""

"""

import numpy as np
                
# compute_velvec2
def compute_velvec2(PLAT,AllVlos,AlldVlos,Allk,AllPlat,AllPlong,Allht,htmin=150.0*1000,htmax=400.0*1000,\
    covar=[1000.*1000.,1000.*1000.,5.*5.],p=[200.0,0.5,2000.0,100.0]):
    
    # 
    # This function takes the line of sight velocities and resolves them into vector velocities as a function of magnetic latitude.  
    # It uses the Bayesian approach and doesn't care about pairs or anything like that.
    # It does care about the magnetic latitude that you choose for binning, so be wise, eh.
    #
    
    # this is the magnetic latitude over which the function bins the data
    plat_in  = np.copy(PLAT)
    Nplout   = plat_in.shape[0]
    plat_out = np.zeros((Nplout),dtype='Float64')
    
    Nparms = Allk.shape[1]
        
    # a priori covariance matrix
    SigmaV = np.matrix(np.diagflat(covar)) 
        
    fracerrs = np.absolute(AlldVlos)/(np.absolute(AllVlos)+p[0])
    abserrs  = np.absolute(AlldVlos)
    
    # loop over output latitudes
    Nmeas    = np.zeros((Nplout))
    chi2     = np.zeros((Nplout))*np.nan
    Vest     = np.zeros((Nplout,Nparms),dtype=AllVlos.dtype)
    dVest    = np.zeros((Nplout,Nparms),dtype=AllVlos.dtype)
    dVestAll = np.zeros((Nplout,Nparms,Nparms),dtype=AllVlos.dtype)
    for i in range(Nplout):
        plat_out[i] = (plat_in[i,0]+plat_in[i,1])/2.0
        I = np.where((AllPlat>plat_in[i,0]) & (AllPlat<plat_in[i,1]) & 
                     (Allht>htmin) & (Allht<htmax) & (np.isfinite(AllVlos)) &
                     (fracerrs<p[1]) & (abserrs<p[3]))[0]
        Vest[i,:]  = np.nan
        dVest[i,:] = np.nan

        if len(I) != 0:
            try:
                tvlos            = np.transpose(np.matrix(AllVlos[I]))
                SigmaE           = np.matrix(np.diagflat(AlldVlos[I]*AlldVlos[I]))
                A                = np.matrix(Allk[I,:])
                tv               = SigmaV*np.transpose(A)*np.linalg.inv(A*SigmaV*np.transpose(A)+SigmaE)*tvlos
                ts               = np.linalg.inv(np.transpose(A)*np.linalg.inv(SigmaE)*A+np.linalg.inv(SigmaV))
                Vest[i,:]        = np.reshape(np.array(tv),(1,Nparms))
                dVest[i,:]       = np.reshape(np.sqrt(np.array(np.diag(ts))),(1,Nparms))
                dVestAll[i,:,:]  = ts
                Nmeas[i]         = len(I)

                temp = np.array(np.matmul(A,Vest[i,:]) - np.transpose(tvlos))
                temp = temp / np.array(AlldVlos[I])

                chi2[i] = np.sum(temp**2) / Nmeas[i]
            except Exception:
                pass
        
    return plat_out, Vest, dVest, dVestAll, Nmeas, chi2
