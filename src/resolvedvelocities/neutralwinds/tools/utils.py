# utils.py

import numpy as np
import warnings

def magnitude_direction(A,Sig,e):
    # Calculate the magnitude of vector A and the clockwise angle between
    # vectors e and A
    # Also calculates corresponding errors
    # A = vector
    # Sig = covariance matrix for A
    # e = vector to take the direction relative to
    # ep = e x z (vector perpendicular to e and up)
    # This is all done with an error analysis using addition in quadrature.
    # All calculations are done with matrix algebra using einsum to prevent
    # nested for loops.
    # Input vectors are assumed to have orthogonal components

    # Calculate the Magnitude of input A

    print(A.shape, Sig.shape, e.shape)

    # Some helper matricies
    # dot product of A and A
    AA = np.einsum('...i,...i->...', A, A)
    # matrix multipy A*Sig*transpose(A)
    ASA = np.einsum('...i,...ij,...j->...', A, Sig, A)

    # calculate magnitude and magnitude error
    magnitude = np.sqrt(AA)
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', category=RuntimeWarning)
        mag_err = np.sqrt(ASA/AA)


    # Now find the angle clockwise around geodetic up in the horizontal plane
    # that is perpendicular to geodetic up, where angle=0 is along the
    # projection 'e' in the horizontal plane

    # Some helper matricies
    # dot product of e and e
    ee = np.einsum('...i,...i->...', e, e)
    # dot product of e and A
    eA = np.einsum('...i,...i->...', e, A)
    # find ep, perpendicular to both e and geodetic up
    ep = np.cross(e,np.array([0,0,1]))
    epep = np.einsum('...i,...i->...', ep, ep)
    epA = np.einsum('...i,...i->...', ep, A)

    # B = ep(e*A)-e(ep*A) = A x (ep x e)
    #B = np.einsum('...ij,...i->...ij',ep,eA)-np.einsum('...ij,...i->...ij',e,epA)
    B = np.einsum('j,...i->...ij',ep,eA)-np.einsum('j,...i->...ij',e,epA)
    # matrix multipy B*Sig*B (covariance propagation)
    BSB = np.einsum('...i,...ij,...j->...', B, Sig, B)

    # calculate direction and direction error
    direction = np.arctan2(np.sqrt(ee)*epA,np.sqrt(epep)*eA)
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', category=RuntimeWarning)
        dir_err = np.sqrt(epep*ee*BSB)/(ee*epA**2-epep*eA**2)

    return magnitude, mag_err, direction*180./np.pi, dir_err*180./np.pi


