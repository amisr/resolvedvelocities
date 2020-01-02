# marp.py

import numpy as np
from apexpy import Apex
import pymap3d as pm

RE = 6371.
hR = 0.

class Marp(Apex):
    def __init__(self, date=None, refh=0, datafile=None, lam0=0., phi0=0., alt=300., geo=False):

        super(Marp, self).__init__(date=date, refh=refh, datafile=None)

        if geo:
            lam0, phi0 = self.geo2apex(lam0, phi0, alt)

        self.lam0 = lam0
        self.phi0 = phi0


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

        # lamr = mlat*np.pi/180.
        # phir = mlon*np.pi/180.
        # lam = alat*np.pi/180.
        # phi = alon*np.pi/180.
        # lam0 = self.lam0*np.pi/180.
        # phi0 = self.phi0*np.pi/180.

        lr = mlat*np.pi/180.
        pr = mlon*np.pi/180.
        l = alat*np.pi/180.
        p = alon*np.pi/180.
        l0 = self.lam0*np.pi/180.
        p0 = self.phi0*np.pi/180.

        f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = self.basevectors_apex(glat, glon, height)
        # print(d1.shape)

        # d1 = np.broadcast_to([1,0,0],(len(lat),3)).T
        # d2 = np.broadcast_to([0,1,0],(len(lat),3)).T
        # d3 = np.broadcast_to([0,0,1],(len(lat),3)).T
        # e1 = np.broadcast_to([1,0,0],(len(lat),3)).T
        # e2 = np.broadcast_to([0,1,0],(len(lat),3)).T
        # e3 = np.broadcast_to([0,0,1],(len(lat),3)).T

        # contravariant derivatives
        dprdp = (np.cos(pr)*np.cos(p-p0) + np.cos(l0)*np.sin(pr)*np.sin(p-p0))*np.cos(l)/np.cos(lr)
        dprdl = (-np.sin(p-p0)*np.cos(pr)*np.sin(l) + np.cos(p-p0)*np.cos(l0)*np.sin(pr)*np.sin(l) - np.sin(l0)*np.sin(pr)*np.cos(l))/np.cos(lr)
        dlrdp = np.sin(l0)*np.sin(p-p0)*np.cos(l)/np.cos(lr)
        dlrdl = (np.sin(l0)*np.sin(l)*np.cos(p-p0) + np.cos(l0)*np.cos(l))/np.cos(lr)
        # dprdp = (np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))*np.cos(l)*np.cos(p-p0)/((np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))**2 + np.sin(p-p0)**2*np.cos(l)**2) + np.sin(p-p0)**2*np.cos(l)**2*np.cos(l0)/((np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))**2 + np.sin(p-p0)**2*np.cos(l)**2)
        # dprdl = -(np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))*np.sin(l)*np.sin(p-p0)/((np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))**2 + np.sin(p-p0)**2*np.cos(l)**2) - (-np.sin(l)*np.cos(l0)*np.cos(p-p0) + np.sin(l0)*np.cos(l))*np.sin(p-p0)*np.cos(l)/((np.sin(l)*np.sin(l0) + np.cos(l)*np.cos(l0)*np.cos(p-p0))**2 + np.sin(p-p0)**2*np.cos(l)**2)
        # dlrdp = np.sin(l0)*np.sin(p-p0)*np.cos(l)/np.sqrt(1 - (np.sin(l)*np.cos(l0) - np.sin(l0)*np.cos(l)*np.cos(p-p0))**2)
        # dlrdl = (np.sin(l)*np.sin(l0)*np.cos(p-p0) + np.cos(l)*np.cos(l0))/np.sqrt(1 - (np.sin(l)*np.cos(l0) - np.sin(l0)*np.cos(l)*np.cos(p-p0))**2)

        # D = 1./np.sqrt(dprdp*dlrdl - dprdl*dlrdp)
        # d1r = (dprdp*d1 + dprdl*d2)*D
        # d2r = (dlrdp*d1 + dlrdl*d2)*D
        # d3r = d3

        sinI = np.abs(2*np.sin(l)/np.sqrt(4-3*np.cos(l)**2))

        # s1 = 1./np.sqrt(dprdp**2/np.cos(l)**2 + dprdl**2/sinI**2)
        # s2 = 1./np.sqrt(dlrdp**2/np.cos(l)**2 + dlrdl**2/sinI**2)

        # d1r = (dprdp*d1/np.cos(l) - dprdl*d2/sinI)/np.sqrt(dprdp**2/np.cos(l)**2 + dprdl**2/sinI**2)
        # d2r = -(dlrdp*d1/np.cos(l) - dlrdl*d2/sinI)/np.sqrt(dlrdp**2/np.cos(l)**2 + dlrdl**2/sinI**2)
        # d3r = d3
        d1r = (dprdp*d1*sinI - dprdl*d2*np.cos(l))/np.sqrt(dprdp**2*sinI**2 + dprdl**2*np.cos(l)**2)
        d2r = -(dlrdp*d1*sinI - dlrdl*d2*np.cos(l))/np.sqrt(dlrdp**2*sinI**2 + dlrdl**2*np.cos(l)**2)
        d3r = d3



        # gradp = d1/((RE+hR)*np.cos(l))
        # gradl = -d2*np.sqrt(4-3*np.cos(l)**2)/(2*(RE+hR)*np.sin(l))
        #
        # gradpr = dprdp*gradp + dprdl*gradl
        # gradlr = dlrdp*gradp + dlrdl*gradl
        #
        # print(d1.shape, gradpr.shape, np.linalg.norm(gradpr, axis=0).shape)
        #
        # d1r = (dprdp*d1 + dprdl*d2)/np.linalg.norm(gradpr, axis=0)
        # d2r = (dlrdp*d1 + dlrdl*d2)/np.linalg.norm(gradlr, axis=0)
        # d3r = d3
        #
        # print(d1r.shape)

        # # covariant derivatives
        # dpdpr = (np.sin(p-p0)*np.cos(l0)*np.sin(pr) + np.cos(p-p0)*np.cos(pr))*np.cos(lr)/np.cos(l)
        # dpdlr = (np.sin(p-p0)*(np.cos(l0)*np.sin(lr)*np.cos(pr) + np.sin(l0)*np.cos(lr)) - np.cos(p-p0)*np.sin(lr)*np.sin(pr))/np.cos(l)
        # dldpr = -np.sin(l0)*np.sin(pr)*np.cos(lr)/np.cos(l)
        # dldlr = (-np.sin(l0)*np.sin(lr)*np.cos(pr) + np.cos(l0)*np.cos(lr))/np.cos(l)




        # dpdpr = (-np.sin(p0)*np.cos(lr)*np.cos(pr) - np.sin(pr)*np.cos(l0)*np.cos(lr)*np.cos(p0))*(np.sin(l0)*np.sin(lr)*np.sin(p0) - np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) - np.sin(pr)*np.cos(lr)*np.cos(p0))/((-np.sin(l0)*np.sin(lr)*np.sin(p0) + np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) + np.sin(pr)*np.cos(lr)*np.cos(p0))**2 + (-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))**2) + (-np.sin(p0)*np.sin(pr)*np.cos(l0)*np.cos(lr) + np.cos(lr)*np.cos(p0)*np.cos(pr))*(-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))/((-np.sin(l0)*np.sin(lr)*np.sin(p0) + np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) + np.sin(pr)*np.cos(lr)*np.cos(p0))**2 + (-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))**2)
        # dpdlr = (np.sin(l0)*np.sin(lr)*np.sin(p0) - np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) - np.sin(pr)*np.cos(lr)*np.cos(p0))*(-np.sin(l0)*np.cos(lr)*np.cos(p0) + np.sin(lr)*np.sin(p0)*np.sin(pr) - np.sin(lr)*np.cos(l0)*np.cos(p0)*np.cos(pr))/((-np.sin(l0)*np.sin(lr)*np.sin(p0) + np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) + np.sin(pr)*np.cos(lr)*np.cos(p0))**2 + (-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))**2) + (-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))*(-np.sin(l0)*np.sin(p0)*np.cos(lr) - np.sin(lr)*np.sin(p0)*np.cos(l0)*np.cos(pr) - np.sin(lr)*np.sin(pr)*np.cos(p0))/((-np.sin(l0)*np.sin(lr)*np.sin(p0) + np.sin(p0)*np.cos(l0)*np.cos(lr)*np.cos(pr) + np.sin(pr)*np.cos(lr)*np.cos(p0))**2 + (-np.sin(l0)*np.sin(lr)*np.cos(p0) - np.sin(p0)*np.sin(pr)*np.cos(lr) + np.cos(l0)*np.cos(lr)*np.cos(p0)*np.cos(pr))**2)
        # dldpr = -np.sin(l0)*np.sin(pr)*np.cos(lr)/np.sqrt(1 - (np.sin(l0)*np.cos(lr)*np.cos(pr) + np.sin(lr)*np.cos(l0))**2)
        # dldlr = (-np.sin(l0)*np.sin(lr)*np.cos(pr) + np.cos(l0)*np.cos(lr))/np.sqrt(1 - (np.sin(l0)*np.cos(lr)*np.cos(pr) + np.sin(lr)*np.cos(l0))**2)

        # D = np.sqrt(dpdpr*dldlr - dldpr*dpdlr)
        # print(D)
        # # D = 1.
        # e1r = (dpdpr*e1 + dldpr*e2)/D
        # e2r = (dpdlr*e1 + dldlr*e2)/D
        # e3r = e3

        e1r = np.cross(d2r.T, d3r.T).T
        e2r = np.cross(d3r.T, d1r.T).T
        e3r = np.cross(d1r.T, d2r.T).T
        #
        # print(e1r.shape)


        return d1r, d2r, d3r, e1r, e2r, e3r
