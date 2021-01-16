# check_reconstruction.py

import numpy as np
import configparser
import pymap3d as pm
from apexpy import Apex
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from resolvedvelocities.synthetic_data.SyntheticData import SyntheticData
from resolvedvelocities.synthetic_data.Field import Field
from resolvedvelocities.synthetic_data.Radar import Radar

def main():

    from argparse import ArgumentParser, RawDescriptionHelpFormatter

    # Build the argument parser tree
    parser = ArgumentParser()
    arg = parser.add_argument('radar',help='Radar name (PFISR or RISRN)')
    arg = parser.add_argument('vvels_config_file',help='Vvels config file.')
    args = vars(parser.parse_args())

    radar_name = args['radar']

    if radar_name=='PFISR':
        # PFISR
        glat = np.arange(62., 71., 1.)
        glon = np.arange(200., 225., 3.)
    elif radar_name=='RISRN':
        # RISRN
        glat = np.arange(73., 83., 1.)
        glon = np.arange(-110., -65., 4.)


    glat, glon = np.meshgrid(glat, glon)
    galt = np.full(glat.shape, 300.)
    glat = glat.flatten()
    glon = glon.flatten()
    galt = galt.flatten()
    field_coords = [glat, glon, galt]

    # generate radar object
    radar = Radar('{}_synth_config.ini'.format(radar_name.lower()))

    trueVe = []
    recstVe = []
    recsterr = []

    rotation_angle = np.arange(0., 360., 15.)

    # for theta in np.arange(0., 2.*np.pi, np.pi/10.):
    for theta in rotation_angle:

        if radar_name=='PFISR':
            # Uniform magnetic field
            Ve1 = 500.*np.cos(theta*np.pi/180.)
            Ve2 = 500.*np.sin(theta*np.pi/180.)
            Ve3 = 0.
            A = Apex(2019)
            f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = A.basevectors_apex(glat, glon, galt)
            field_values = Ve1*e1.T + Ve2*e2.T + Ve3*e3.T

        elif radar_name=='RISRN':
            # Uniform field
            e = 500.*np.cos(theta*np.pi/180.)
            n = 500.*np.sin(theta*np.pi/180.)
            u = 0.
            vx, vy, vz = pm.enu2uvw(e, n, u, np.mean(glat), np.mean(glon))
            ve, vn, vu = pm.uvw2enu(vx, vy, vz, glat, glon)
            field_values = np.array([ve, vn, vu]).T


        # generate field object
        field = Field(apex_year=2019, field_coords=field_coords, field_values=field_values)

        # use field and radar objects to produce synthetic data set
        synth_data = SyntheticData(field, radar)

        # run resolvevectors algothrithm on synthetic data set
        rv = synth_data.eval_vvels(args['vvels_config_file'])


        targalt = 300.
        aidx = np.argmin(np.abs(rv.bin_galt[:,0]-targalt))

        # convert bin locations to ECEF
        X, Y, Z = pm.geodetic2ecef(rv.bin_glat[aidx,:], rv.bin_glon[aidx,:], rv.bin_galt[aidx,:]*1000.)
        # interpolate field to bin locations
        Vx = field.interpVx(np.array([X,Y,Z]).T)
        Vy = field.interpVy(np.array([X,Y,Z]).T)
        Vz = field.interpVz(np.array([X,Y,Z]).T)
        # convert field components to marp
        Ve, Vn, Vu = pm.uvw2enu(Vx, Vy, Vz, rv.bin_glat[aidx,:], rv.bin_glon[aidx,:])
        d1,d2,d3,e1,e2,e3 = rv.marp.basevectors_marp(rv.bin_glat[aidx,:], rv.bin_glon[aidx,:], rv.bin_galt[aidx,:])
        Ve1 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d1.T)
        Ve2 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d2.T)
        Ve3 = np.einsum('...i,...i->...',np.array([Ve, Vn, Vu]).T,d3.T)


        trueVe.append(np.array([Ve1,Ve2,Ve3]).T)
        recstVe.append(rv.Velocity[0,:,:])
        recsterr.append(np.array([np.sqrt(rv.VelocityCovariance[0,:,i,i]) for i in range(3)]).T)


    trueVe = np.array(trueVe)
    recstVe = np.array(recstVe)
    recsterr = np.array(recsterr)

    # form array of bin numbers
    bins = np.arange(len(rv.bin_glat[aidx,:]))+1

    # create summary plot
    fig = plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(3,4)
    gs.update(left=0.05,right=0.95,bottom=0.05,top=0.95,hspace=0.3)

    for i in range(3):
        ax = plt.subplot(gs[i,0])
        c = ax.pcolormesh(rotation_angle, bins, trueVe[:,:,i].T, vmin=-500., vmax=500., cmap=plt.get_cmap('rainbow'))
        ax.set_title('True Ve{}'.format(i+1))
        ax.set_xlabel('Rotation Angle ($^\circ$)')
        ax.set_ylabel('Bin #')
        plt.colorbar(c)

        ax = plt.subplot(gs[i,1])
        c = ax.pcolormesh(rotation_angle, bins, recstVe[:,:,i].T, vmin=-500., vmax=500., cmap=plt.get_cmap('rainbow'))
        ax.set_title('Reconstructed Ve{}'.format(i+1))
        ax.set_xlabel('Rotation Angle ($^\circ$)')
        ax.set_ylabel('Bin #')
        plt.colorbar(c)

        ax = plt.subplot(gs[i,2])
        c = ax.pcolormesh(rotation_angle, bins, recsterr[:,:,i].T, vmin=0., vmax=100., cmap=plt.get_cmap('Reds'))
        ax.set_title('Reconstructed Ve{} Error'.format(i+1))
        ax.set_xlabel('Rotation Angle ($^\circ$)')
        ax.set_ylabel('Bin #')
        plt.colorbar(c)

        ax = plt.subplot(gs[i,3])
        c = ax.pcolormesh(rotation_angle, bins, np.abs(trueVe[:,:,i].T-recstVe[:,:,i].T), vmin=0., vmax=100., cmap=plt.get_cmap('Reds'))
        ax.set_title('|True Ve{0} - Reconstructed Ve{0}|'.format(i+1))
        ax.set_xlabel('Rotation Angle ($^\circ$)')
        ax.set_ylabel('Bin #')
        plt.colorbar(c)

    plt.show()




if __name__ == '__main__':
    main()
