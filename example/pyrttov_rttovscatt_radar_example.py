
# Example program calling RTTOV-SCATT radar simulation using pyrttov
# Requires the file lib/rttov_wrapper_f2py.so to be in your $PYTHONPATH or the current directory

# Example profile data is contained in example_data_rttovscatt.py

import pyrttov
from example import example_data_rttovscatt as ex

import numpy as np
import sys
# import matplotlib
import pylab as plt

rttov_installdir = '/home/htht/rttov/rttov13/'

if __name__ == '__main__':

    # Run direct (0) or K (1)
    doK=1

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(ex.p_ex)
    nprofiles = 2
    myProfiles = pyrttov.ProfilesScatt(nprofiles, nlevels)

    # Associate the profiles and other data from example_data_rttovscatt.h with myProfiles
    # Note that the zeeman data are not mandatory and are omitted here

    def expand2nprofiles(n, nprof):
        """Transform 1D array to a [nprof, nlevels] array"""
        outp = np.empty((nprof, len(n)), dtype=n.dtype)
        for i in range(nprof):
            outp[i, :] = n[:]
        return outp

    myProfiles.GasUnits = ex.gas_units
    myProfiles.P = expand2nprofiles(ex.p_ex, nprofiles)
    myProfiles.T = expand2nprofiles(ex.t_ex, nprofiles)
    myProfiles.Q = expand2nprofiles(ex.q_ex, nprofiles)

    #myProfiles.Ph = expand2nprofiles(ex.ph_ex, nprofiles)
    ph = np.zeros((nprofiles,nlevels+1), dtype=np.float64)
    for i in range(nprofiles):
        ph[i,:nlevels] = ex.ph_ex
        ph[i,nlevels] = ex.s2m[i,0]
    myProfiles.Ph = ph

    # Modified cloud cover to represent a rain fraction below the peak in the rain-producing cloud
    modified_cc = ex.cc_ex
    modified_cc[53:] = modified_cc[52]
    myProfiles.HydroFrac = expand2nprofiles(modified_cc, nprofiles)
    #myProfiles.HydroFrac = expand2nprofiles(ex.cc_ex, nprofiles)
    #myProfiles.HydroFrac = np.ones((nprofiles, len(ex.cc_ex)), dtype=ex.cc_ex.dtype) # AJGDB to replace some extensive changes in rttov_hydro.F90

    myProfiles.Clw = expand2nprofiles(ex.clw_ex, nprofiles)
    myProfiles.Ciw = expand2nprofiles(ex.ciw_ex, nprofiles)
    myProfiles.Snow = expand2nprofiles(ex.snow_ex, nprofiles)
    myProfiles.Rain = expand2nprofiles(ex.rain_ex, nprofiles)

    myProfiles.Angles = ex.angles
    myProfiles.S2m = ex.s2m
    myProfiles.Skin = ex.skin
    myProfiles.SurfType = ex.surftype
    myProfiles.SurfGeom = ex.surfgeom
    myProfiles.DateTimes = ex.datetimes

    # ------------------------------------------------------------------------
    # Set up RttovScatt instance
    # ------------------------------------------------------------------------

    dprRttov = pyrttov.RttovScatt()

    dprRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/rttov13pred54L/rtcoef_gpm_1_dpr.dat")
    dprRttov.FileHydrotable = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/hydrotable/hydrotable_gpm_dpr.dat")

    dprRttov.Options.LuserCfrac = False
    dprRttov.Options.VerboseWrapper = True
    dprRttov.Options.StoreRad = True

    #Activate reflectivity calculations
    dprRttov.CalcZef = True

    # Load instrument
    try:
        dprRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument: {!s}\n".format(e))
        sys.exit(1)

    # Associate the profiles with RttovScatt instance
    try:
        dprRttov.Profiles = myProfiles
    except pyrttov.RttovError as e:
        sys.stderr.write("Error setting profiles: {!s}\n".format(e))
        sys.exit(1)

    nchan_dpr = dprRttov.Nchannels

    # Set up the surface emissivity array and associate with the RttovScatt object
    surfemis_dpr = np.zeros((nprofiles,nchan_dpr), dtype=np.float64)

    surfemis_dpr[:,:] = -1.
    dprRttov.SurfEmis = surfemis_dpr

    # ------------------------------------------------------------------------
    # Call RTTOV-SCATT
    # ------------------------------------------------------------------------


    try:
        if doK:
            dprRttov.runK()
        else:
            dprRttov.runDirect()

    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV-SCATT: {!s}\n".format(e))
        sys.exit(1)

    # ------------------------------------------------------------------------
    # Plot some of the output
    # ------------------------------------------------------------------------

    print(dprRttov.Zef[0, 0, :])
    print(' ')
    print(dprRttov.AZef[0, 0, :])
    print(' ')

    # Use the altitude of Ku from the first profile (since it is identical between channels and profiles anyway)
    altitude_km = dprRttov.GeometricHeight[0, 0, :]/1000

    if 1:
        plt.plot(dprRttov.Zef[0, 0, :], altitude_km, 'b', label='Ku')
        plt.plot(dprRttov.AZef[0, 0, :], altitude_km, '--b', label='Ku Attenuated')
        plt.plot(dprRttov.Zef[0, 1, :], altitude_km, 'r', label='Ka')
        plt.plot(dprRttov.AZef[0, 1, :], altitude_km, '--r', label='Ka Attenuated')
        plt.legend()
        plt.xlim([-15, 25])
        plt.ylim([0., 8.])
        plt.xlabel('Reflectitivy (dBz)')
        plt.ylabel('Height (km)')
        plt.show()

    if doK:
        plt.plot(dprRttov.GasesK[2, 0, 0, :], altitude_km, 'y', label='CLW')
        plt.plot(dprRttov.GasesK[3, 0, 0, :], altitude_km, 'g', label='CIW')
        plt.plot(dprRttov.GasesK[4, 0, 0, :], altitude_km, 'b', label='RAIN')
        plt.plot(dprRttov.GasesK[5, 0, 0, :], altitude_km, 'r', label='SNOW')
        plt.legend()
        plt.ylim([0., 12.])
        plt.xlabel('Hydrometeor Jacobians (dBz/(kg/kg))')
        plt.ylabel('Height (km)')
        plt.show()

    plt.plot(myProfiles.Clw[0,:],  altitude_km, 'y', label='CLW')
    plt.plot(myProfiles.Ciw[0,:],  altitude_km, 'g', label='CIW')
    plt.plot(myProfiles.Rain[0,:], altitude_km, 'b', label='RAIN')
    plt.plot(myProfiles.Snow[0,:],   altitude_km, 'r', label='SNOW')
    plt.legend()
    plt.ylim([0., 12.])
    plt.xlabel('Hydrometeor profiles ((kg/kg))')
    plt.ylabel('Height (km)')
    plt.show()

    plt.plot(myProfiles.HydroFrac[0,:],  altitude_km, 'y', label='CC')
    plt.legend()
    plt.ylim([0., 12.])
    plt.xlabel('Hydrometeor profiles ((kg/kg))')
    plt.ylabel('Height (km)')
    plt.show()
