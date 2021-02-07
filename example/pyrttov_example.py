
# Example of using the Rttov class to call RTTOV for multiple instruments
# with the emissivity and BRDF atlases

# Three Rttov instances are created representing three instruments

import pyrttov
from example import example_data as ex
import numpy as np
import sys

rttov_installdir = r'/home/htht/rttov/rttov13/'

if __name__ == '__main__':

    # This example program simulates two profiles for each of three instruments
    # The example profile data are defined in example_data

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(ex.p_ex)

    nprofiles = 2
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # Associate the profiles and other data from example_data.h with myProfiles
    # Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and
    # are omitted here

    def expand2nprofiles(n, nprof):
        """Transform 1D array to a [nprof, nlevels] array"""
        outp = np.empty((nprof, len(n)), dtype=n.dtype)
        for i in range(nprof):
            outp[i, :] = n[:]
        return outp

    myProfiles.GasUnits = ex.gas_units
    myProfiles.P = expand2nprofiles(ex.p_ex, nprofiles)
    myProfiles.T = expand2nprofiles(ex.t_ex, nprofiles)
    # Modify the temperature of the second profile
    myProfiles.T[1, :] += 2
    myProfiles.Q = expand2nprofiles(ex.q_ex, nprofiles)
    myProfiles.CO2 = expand2nprofiles(ex.co2_ex, nprofiles)

    myProfiles.Angles = ex.angles
    myProfiles.S2m = ex.s2m
    myProfiles.Skin = ex.skin
    myProfiles.SurfType = ex.surftype
    myProfiles.SurfGeom = ex.surfgeom
    myProfiles.DateTimes = ex.datetimes

    # ------------------------------------------------------------------------
    # Set up Rttov instances for each instrument
    # ------------------------------------------------------------------------

    # Create three Rttov objects for three instruments
    seviriRttov = pyrttov.Rttov()
    hirsRttov = pyrttov.Rttov()
    mhsRttov = pyrttov.Rttov()

    # For HIRS and MHS we will read all channels, but we will read a subset
    # for SEVIRI
    nchan_seviri = 10
    nchan_hirs = 19
    nchan_mhs = 5

    # For SEVIRI exclude ozone and hi-res vis channels (9 and 12) in this
    # example
    chan_list_seviri = (1, 2, 3, 4, 5, 6, 7, 9, 10, 11)

    # Set the options for each Rttov instance:
    # - the path to the coefficient file must always be specified
    # - turn RTTOV interpolation on (because input pressure levels differ from
    #   coefficient file levels)
    # - set the verbose_wrapper flag to true so the wrapper provides more
    #   information
    # - enable solar simulations for SEVIRI
    # - enable CO2 simulations for HIRS (the CO2 profiles are ignored for
    #   the SEVIRI and MHS simulations)
    # - enable the store_trans wrapper option for MHS to provide access to
    #   RTTOV transmission structure

    seviriRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/rttov13pred54L/rtcoef_msg_3_seviri_o3.dat")
    seviriRttov.Options.AddInterp = True
    seviriRttov.Options.AddSolar = True
    seviriRttov.Options.VerboseWrapper = True

    hirsRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                        "rtcoef_rttov13/rttov13pred54L/rtcoef_metop_1_hirs_o3co2.dat")
    hirsRttov.Options.AddInterp = True
    hirsRttov.Options.CO2Data = True
    hirsRttov.Options.VerboseWrapper = True

    mhsRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                       "rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_19_mhs.dat")
    mhsRttov.Options.AddInterp = True
    mhsRttov.Options.StoreTrans = True
    mhsRttov.Options.VerboseWrapper = True

    # Load the instruments: for HIRS and MHS do not supply a channel list and
    # so read all channels
    try:
        seviriRttov.loadInst(chan_list_seviri)
        hirsRttov.loadInst()
        mhsRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    # Associate the profiles with each Rttov instance
    seviriRttov.Profiles = myProfiles
    hirsRttov.Profiles = myProfiles
    mhsRttov.Profiles = myProfiles

    # ------------------------------------------------------------------------
    # Load the emissivity and BRDF atlases
    # ------------------------------------------------------------------------

    # Load the emissivity and BRDF atlases:
    # - load data for the month in the profile data
    # - load the IR emissivity atlas data for multiple instruments so it can be used for SEVIRI and HIRS
    # - SEVIRI is the only VIS/NIR instrument we can use the single-instrument initialisation for the BRDF atlas

    irAtlas = pyrttov.Atlas()
    irAtlas.AtlasPath = r'/home/htht/rttov/rttov13/emis_data/IR/'
    irAtlas.loadIrEmisAtlas(ex.datetimes[0][1], ang_corr=True) # Include angular correction, but do not initialise for single-instrument

    brdfAtlas = pyrttov.Atlas()
    brdfAtlas.AtlasPath = r'/home/htht/rttov/rttov13/brdf_data/'
    brdfAtlas.loadBrdfAtlas(ex.datetimes[0][1], seviriRttov) # Supply Rttov object to enable single-instrument initialisation
    brdfAtlas.IncSea = False                                 # Do not use BRDF atlas for sea surface types

    # TELSEM2 atlas does not require an Rttov object to initialise
    mwAtlas = pyrttov.Atlas()
    mwAtlas.AtlasPath = r'/home/htht/rttov/rttov13/emis_data/MW/'
    mwAtlas.loadMwEmisAtlas(ex.datetimes[0][1])

    # Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    surfemisrefl_seviri = np.zeros((4,nprofiles,nchan_seviri), dtype=np.float64)
    surfemisrefl_hirs = np.zeros((4,nprofiles,nchan_hirs), dtype=np.float64)
    surfemisrefl_mhs = np.zeros((4,nprofiles,nchan_mhs), dtype=np.float64)

    seviriRttov.SurfEmisRefl = surfemisrefl_seviri
    hirsRttov.SurfEmisRefl = surfemisrefl_hirs
    mhsRttov.SurfEmisRefl = surfemisrefl_mhs

    # ------------------------------------------------------------------------
    # Call RTTOV
    # ------------------------------------------------------------------------

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
    # calcemis/calcrefl TRUE - see RTTOV user guide)

    surfemisrefl_seviri[:,:,:] = -1.
    surfemisrefl_hirs[:,:,:]   = -1.
    surfemisrefl_mhs[:,:,:]    = -1.

    # Call emissivity and BRDF atlases
    try:
        # Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
        # *loaded* channels which is what is required
        surfemisrefl_seviri[0,:,:] = irAtlas.getEmisBrdf(seviriRttov)
        surfemisrefl_seviri[1,:,:] = brdfAtlas.getEmisBrdf(seviriRttov)
        surfemisrefl_hirs[0,:,:] = irAtlas.getEmisBrdf(hirsRttov)
        surfemisrefl_mhs[0,:,:] = mwAtlas.getEmisBrdf(mhsRttov)

    except pyrttov.RttovError as e:
        # If there was an error the emissivities/BRDFs will not have been modified so it
        # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        sys.stderr.write("Error calling atlas: {!s}".format(e))

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are
    # simulated
    try:
        seviriRttov.runDirect()
        hirsRttov.runDirect()
        mhsRttov.runDirect()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)

    # ------------------------------------------------------------------------
    # Print out some of the output
    # ------------------------------------------------------------------------

    print
    print("SELECTED OUTPUT")
    print

    print("SEVIRI visible channel reflectances, channels 1-3")
    for p in range(nprofiles):
        print("Profile {:d}:".format(p))
        for c in range(3):
            print("  Ch #{:02d} refl={:f}".format(chan_list_seviri[c],
                                                  seviriRttov.BtRefl[p, c]))
        print

    print("HIRS radiances")
    for p in range(nprofiles):
        print("Profile {:d}:".format(p))
        for c in range(nchan_hirs):
            print("  Ch #{:02d} rad={:f}".format(c + 1, hirsRttov.Rads[p, c]))
        print

    # We can access the RTTOV transmission structure because the store_trans
    # option was set above for mhsRttov
    print("MHS total transmittance")
    for p in range(nprofiles):
        print("Profile {:d}:".format(p))
        for c in range(nchan_mhs):
            print("  Ch #{:02d} tau={:f}".format(c + 1,
                                                 mhsRttov.TauTotal[p, c]))
        print

    # ------------------------------------------------------------------------
    # Deallocate memory
    # ------------------------------------------------------------------------


    # Because of Python's garbage collector, there should be no need to
    # explicitly deallocate memory
