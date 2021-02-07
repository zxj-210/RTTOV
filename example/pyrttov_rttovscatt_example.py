
# Example program calling RTTOV-SCATT using pyrttov
# Requires the file lib/rttov_wrapper_f2py.so to be in your $PYTHONPATH or the current directory

# Example profile data is contained in example_data_rttovscatt.py

import pyrttov
from example import example_data_rttovscatt as ex

import numpy as np
import sys

rttov_installdir = '../'

if __name__ == '__main__':

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

    myProfiles.HydroFrac = expand2nprofiles(ex.cc_ex, nprofiles)
    myProfiles.Clw = expand2nprofiles(ex.clw_ex, nprofiles)
    myProfiles.Ciw = expand2nprofiles(ex.ciw_ex, nprofiles)
    myProfiles.Snow = expand2nprofiles(ex.snow_ex, nprofiles)
    myProfiles.Rain = expand2nprofiles(ex.rain_ex, nprofiles)

    # Equivalent code using the flexible hydro interface
    #myProfiles.HydroFrac1 = expand2nprofiles(ex.cc_ex, nprofiles)
    #myProfiles.Hydro4 = expand2nprofiles(ex.clw_ex, nprofiles)
    #myProfiles.Hydro5 = expand2nprofiles(ex.ciw_ex, nprofiles)
    #myProfiles.Hydro2 = expand2nprofiles(ex.snow_ex, nprofiles)
    #myProfiles.Hydro1 = expand2nprofiles(ex.rain_ex, nprofiles)

    myProfiles.Angles = ex.angles
    myProfiles.S2m = ex.s2m
    myProfiles.Skin = ex.skin
    myProfiles.SurfType = ex.surftype
    myProfiles.SurfGeom = ex.surfgeom
    myProfiles.DateTimes = ex.datetimes


    # ------------------------------------------------------------------------
    # Set up RttovScatt instance
    # ------------------------------------------------------------------------

    amsuaRttov = pyrttov.RttovScatt()

    amsuaRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_15_amsua.dat")
    amsuaRttov.FileHydrotable = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/hydrotable/hydrotable_noaa_amsua.dat")

    amsuaRttov.Options.LuserCfrac = False
    amsuaRttov.Options.VerboseWrapper = True
    amsuaRttov.Options.StoreRad = True

    # Load instrument
    try:
        amsuaRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument: {!s}\n".format(e))
        sys.exit(1)

    # Associate the profiles with RttovScatt instance
    try:
        amsuaRttov.Profiles = myProfiles
    except pyrttov.RttovError as e:
        sys.stderr.write("Error setting profiles: {!s}\n".format(e))
        sys.exit(1)

    nchan_amsua = amsuaRttov.Nchannels


    # ------------------------------------------------------------------------
    # Load the emissivity atlas
    # ------------------------------------------------------------------------

    # TELSEM2 atlas does not require an Rttov object to initialise
    mwAtlas = pyrttov.Atlas()
    mwAtlas.AtlasPath = '{}/{}'.format(rttov_installdir, "emis_data")
    mwAtlas.loadMwEmisAtlas(ex.datetimes[0][1])
    mwAtlas.IncSea = False

    # Set up the surface emissivity array and associate with the RttovScatt object
    surfemis_amsua = np.zeros((nprofiles,nchan_amsua), dtype=np.float64)

    amsuaRttov.SurfEmis = surfemis_amsua

    # ------------------------------------------------------------------------
    # Call RTTOV-SCATT
    # ------------------------------------------------------------------------

    # Surface emissivity array must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity values (i.e. equivalent to
    # calcemis TRUE - see RTTOV user guide)

    surfemis_amsua[:,:] = -1.

    # Call emissivity atlas
    try:
        surfemis_amsua[:,:] = mwAtlas.getEmisBrdf(amsuaRttov)
    except pyrttov.RttovError as e:
        # If there was an error the emissivities/BRDFs will not have been modified so it
        # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        sys.stderr.write("Error calling atlas: {!s}\n".format(e))

    try:
        #amsuaRttov.runDirect()
        amsuaRttov.runK()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}\n".format(e))
        sys.exit(1)


    # ------------------------------------------------------------------------
    # Print out some of the output
    # ------------------------------------------------------------------------

    print('Surface emissivity used by RTTOV')
    print(amsuaRttov.SurfEmis[:,:])

    print('Total cloudy BT')
    print(amsuaRttov.Bt)

    print('Clear-sky BT')
    print(amsuaRttov.BtClear)
    
    # Example code for rain Jacobian using flexible hydro interface
    #   (profiles must have been specified using the same interface)
    #print('Rain Jacobian for profile 1, channel 1')
    #print(amsuaRttov.RainK[0,0,:])
    #print(amsuaRttov.getHydroNK(1)[0,0,:])

