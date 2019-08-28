# marp.py

import numpy as np
from apexpy import Apex
import pymap3d as pm

RE = 6371.
hR = 0.

class Marp(Apex):
    def __init__(self, date=None, refh=0, datafile=None, lam0=0., phi0=0.):

        self.lam0 = lam0
        self.phi0 = phi0

        super(Marp, self).__init__(date=date, refh=refh, datafile=None)

    def apex2marp(self, lam, phi):

        lam = lam*np.pi/180.
        phi = phi*np.pi/180.
        lam0 = self.lam0*np.pi/180.
        phi0 = self.phi0*np.pi/180.

        xr = np.cos(lam0)*np.cos(lam)*np.cos(phi-phi0) + np.sin(lam0)*np.sin(lam)
        yr = np.cos(lam)*np.sin(phi-phi0)
        zr = -np.sin(lam0)*np.cos(lam)*np.cos(phi-phi0) + np.cos(lam0)*np.sin(lam)

        phir = np.arctan2(yr, xr)
        lamr = np.arcsin(zr)

        return lamr*180./np.pi, phir*180./np.pi

    def marp2apex(self, lamr, phir):

        lamr = lamr*np.pi/180.
        phir = phir*np.pi/180.
        lam0 = self.lam0*np.pi/180.
        phi0 = self.phi0*np.pi/180.

        x = np.cos(phi0)*np.cos(lam0)*np.cos(lamr)*np.cos(phir) - np.cos(phi0)*np.sin(lam0)*np.sin(lamr) - np.sin(phi0)*np.cos(lamr)*np.sin(phir)
        y = np.sin(phi0)*np.cos(lam0)*np.cos(lamr)*np.cos(phir) - np.sin(phi0)*np.sin(lam0)*np.sin(lamr) + np.cos(phi0)*np.cos(lamr)*np.sin(phir)
        z = np.sin(lam0)*np.cos(lamr)*np.cos(phir) + np.cos(lam0)*np.sin(lamr)

        phi = np.arctan2(y, x)
        lam = np.arcsin(z)

        return lam*180./np.pi, phi*180./np.pi

    def geo2marp(self, glat, glon, height):
        alat, alon = self.geo2apex(glat, glon, height)
        mlat, mlon = self.apex2marp(alat, alon)
        return mlat, mlon

    def marp2geo(self, mlat, mlon, height):
        alat, alon = self.marp2apex(mlat, mlon)
        glat, glon, err = self.apex2geo(alat, alon, height)
        return glat, glon, err

    def basevectors_marp(self, lat, lon, height, coords='geo'):

        if coords == 'geo':
            glat = lat
            glon = lon
            alat, alon = self.geo2apex(glat, glon, height)
            mlat, mlon = self.apex2marp(alat, alon)

        if coords == 'apex':
            alat = lat
            alon = lon
            glat, glon, _ = self.apex2geo(alat, alon, height)
            mlat, mlon = self.apex2marp(alat, alon)

        lamr = mlat*np.pi/180.
        phir = mlon*np.pi/180.
        lam = alat*np.pi/180.
        phi = alon*np.pi/180.
        lam0 = self.lam0*np.pi/180.
        phi0 = self.phi0*np.pi/180.

        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.basevectors_apex(glat, glon, height)

        dprdp = (np.cos(phir)*np.cos(phi-phi0) + np.cos(lam0)*np.sin(phir)*np.sin(phi-phi0))*np.cos(lam)/np.cos(lamr)
        dprdl = (-np.sin(phi-phi0)*np.cos(phir)*np.sin(lam) + np.cos(phi-phi0)*np.cos(lam0)*np.sin(phir)*np.sin(lam) - np.sin(lam0)*np.sin(phir)*np.cos(lam))/np.cos(lamr)
        dlrdp = np.sin(lam0)*np.sin(phi-phi0)*np.cos(lam)/np.cos(lamr)
        dlrdl = (np.sin(lam0)*np.sin(lam)*np.cos(phi-phi0) + np.cos(lam0)*np.cos(lam))/np.cos(lamr)
        D = 1./np.sqrt(dprdp*dlrdl - dprdl*dlrdp)

        d1r = (dprdp*d1 + dprdl*d2)*D
        d2r = (dlrdp*d1 + dlrdl*d2)*D
        d3r = d3


        dpdpr = (np.sin(phi-phi0)*np.cos(lam0)*np.sin(phir) + np.cos(phi-phi0)*np.cos(phir))*np.cos(lamr)/np.cos(lam)
        dpdlr = (np.sin(phi-phi0)*(np.cos(lam0)*np.sin(lamr)*np.cos(phir) + np.sin(lam0)*np.cos(lamr)) - np.cos(phi-phi0)*np.sin(lamr)*np.sin(phir))/np.cos(lam)
        dldpr = -np.sin(lam0)*np.sin(phir)*np.cos(lamr)/np.cos(lam)
        dldlr = (-np.sin(lam0)*np.sin(lamr)*np.cos(phir) + np.cos(lam0)*np.cos(lamr))/np.cos(lam)
        D = np.sqrt(dpdpr*dldlr - dldpr*dpdlr)

        e1r = (dpdpr*e1 + dldpr*e2)/D
        e2r = (dpdlr*e1 + dldlr*e2)/D
        e3r = e3

        return d1r, d2r, d3r, e1r, e2r, e3r

