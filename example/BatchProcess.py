#coding:utf-8

import os
import sys
import datetime
import glob
#import h5py
import numpy as np
from ncpro import lb_writenc, lb_readnc

if __name__ == '__main__':

    pathin = r'/home/htht/rttov/prog/data'
    pathout = r'/home/htht/rttov/prog/result'

    albedo_vis = []
    albedo_short = []
    albedo_nir = []

    filelist = glob.glob(os.path.join(pathin, '*20210115*'))

    for filename in filelist :
        name = os.path.basename(filename)

        # cmd = '{exe} {inname} {outname}'.format(
        #     exe = '/home/htht/rttov/prog/albedo/albedo.e',
        #     inname = filename,
        #     outname = os.path.join(pathout, name.replace('_L1_GEO-_','_L2_LSA-_'))
        # )
        #
        # print(cmd)
        # os.system(cmd)

        outname = os.path.join(pathout, name.replace('_L1-_GEO-_','_L2-_LSA-_'))

        #with h5py.File(outname, 'r') as fp :
        albedo_vis.append(lb_readnc(outname, 'albedo_vis'))
        albedo_nir.append(lb_readnc(outname, 'albedo_nir'))
        albedo_short.append(lb_readnc(outname, 'albedo_short'))

    albedo_vis = np.array(albedo_vis)
    albedo_short = np.array(albedo_short)
    albedo_nir = np.array(albedo_nir)

    albedo_vis[albedo_vis > 1 or albedo_vis < 0] = np.nan
    albedo_nir[albedo_nir > 1 or albedo_nir < 0] = np.nan
    albedo_short[albedo_short > 1 or albedo_short < 0] = np.nan

    albedo_vis = np.nanmax(albedo_vis, axis=0)
    albedo_nir = np.nanmax(albedo_nir, axis=0)
    albedo_short = np.nanmax(albedo_short, axis=0)

    albedo_vis[np.isnan(albedo_vis)] = 65535.0
    albedo_nir[np.isnan(albedo_nir)] = 65535.0
    albedo_short[np.isnan(albedo_short)] = 65535.0

    outname = 'FY4A-_AGRI--_N_DISK_1047E_L2-_LSA-_MULT_NOM_20210115000000_20210115235959_4000M_V0001.NC'

    lb_writenc(outname, 'albedo_vis', albedo_vis, dimension=('y', 'x'), overwrite=1)
    lb_writenc(outname, 'albedo_nir', albedo_nir, dimension=('y', 'x'), overwrite=0)
    lb_writenc(outname, 'albedo_short', albedo_short, dimension=('y', 'x'), overwrite=0)