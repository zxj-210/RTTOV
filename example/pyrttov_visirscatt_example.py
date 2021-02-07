
# Example of using the Rttov class to do visible/IR scattering simulations
# with explicit optical parameters

import pyrttov
from example import example_data as ex, example_data_opt_param as exop

import numpy as np
import sys

rttov_installdir = '../'

if __name__ == '__main__':

    # This example program simulates two profiles
    # The example profile data are defined in example_data
    # The example optical properties are defined in example_data_opt_param

    # ------------------------------------------------------------------------
    # Set up the profile data
    # ------------------------------------------------------------------------

    # Declare an instance of Profiles
    nlevels = len(ex.p_ex)
    nprofiles = 2
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # Associate the profiles and other data from example_data with myProfiles
    # (The optical parameter data is handled below)

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

    myProfiles.Angles = ex.angles
    myProfiles.S2m = ex.s2m
    myProfiles.Skin = ex.skin
    myProfiles.SurfType = ex.surftype
    myProfiles.SurfGeom = ex.surfgeom
    myProfiles.DateTimes = ex.datetimes

    # ------------------------------------------------------------------------
    # Set up the Rttov instance
    # ------------------------------------------------------------------------

    seviriRttov = pyrttov.Rttov()

    # Set the options
    seviriRttov.FileCoef = '{}/{}'.format(rttov_installdir,
                                          "rtcoef_rttov13/rttov13pred54L/rtcoef_msg_3_seviri_o3.dat")
    seviriRttov.Options.AddInterp = True        # Use the RTTOV interpolator
    seviriRttov.Options.AddSolar = True         # Turn on solar radiation
    seviriRttov.Options.VerboseWrapper = True   # Turn on verbose wrapper output
    seviriRttov.Options.AddClouds = True        # Activate cloud simulations
    seviriRttov.Options.UserCldOptParam = True  # Use explicit cloud optical properties
    seviriRttov.Options.VisScattModel = 1       # Use DOM solver for solar radiation
    seviriRttov.Options.IrScattModel = 2        # Use Chou-scaling for thermal emitted radiation
    seviriRttov.Options.DomNstreams = 16        # Number of DOM streams to use
    seviriRttov.Options.Nthreads = 8            # Take advantage of multiple threads if RTTOV was compiled with OpenMP
    seviriRttov.Options.StoreRad = True         # Store all radiance outputs

    # Load the instrument (reads all channels)
    try:
        seviriRttov.loadInst()
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument: {!s}".format(e))
        sys.exit(1)

    # ------------------------------------------------------------------------
    # Specify cloud optical parameters
    # ------------------------------------------------------------------------

    # The example_data_opt_param file contains data for a single cloud layer for two channels:
    #   opt_param_chan_list[nchan]  - the channels for which optical parameters are provided
    #   abscoef[nchan]              - absorption coefficient for each channel
    #   scacoef[nchan]              - scattering coefficient for each channel
    #   phangle[nphangle]           - phase function angle grid
    #   phasefn[nchan][nphangle]    - phase functions for each channel

    # In this example we define a single-layer cloud for two channels: one visible, one IR
    chan_list = exop.opt_param_chan_list
    nchan = len(chan_list)

    # For cloudy simulations the cloud fraction profile must be specified in the Profiles object
    cfrac = np.zeros((nprofiles, nlevels), dtype=np.float64)
    cfrac[:,74] = 1.                      # cloud fraction: both profiles, layer 75
    myProfiles.Cfrac = cfrac

    # Define optical parameters:
    # - these are defined for every *layer* for every channel being simulated for every profile
    # - absorption and scattering coefficients are always required
    # - bpr parameter is only required when using Chou-scaling, can be zero otherwise; this can
    #   be calculated from the phase function using the calcBpr method of Rttov
    # - phase functions are required only for solar-affected channels when the AddSolar option is true
    # - Legendre coefficients are only required if using the DOM solver; these can be calculated
    #   from the phase function using the calcLegcoef method of Rttov

    # In the second profile the scattering coefficient is doubled

    # Absorption coef, scattering coef, bpr parameter [3][nprofiles][nchan][nlayers]
    asb = np.zeros((3,nprofiles,nchan,nlevels-1), dtype=np.float64)
    asb[0,0,:,74] = exop.abscoef[:]       # abs coef: profile 1, both channels, layer 75
    asb[0,1,:,74] = exop.abscoef[:]       # abs coef: profile 2, both channels, layer 75
    asb[1,0,:,74] = exop.scacoef[:]       # sca coef: profile 1, both channels, layer 75
    asb[1,1,:,74] = exop.scacoef[:] * 2.  # sca coef: profile 2, both channels, layer 75

    # Since Chou-scaling is used in the IR we must calculate the bpr for the IR channel.
    # Note that this is relatively slow so for performance-critical applications this
    # calculation would be done off-line for each phase function. If RTTOV was compiled
    # with OpenMP it uses the number of threads specified in the wrapper options.
    # The bpr value is left as zero for the visible channel.
    bpr = seviriRttov.calcBpr(exop.phangle, exop.phasefn[1,:])
    asb[2,:,1,74] = bpr                   # bpr: both profiles, channel 2, layer 75

    # Specify phase function for the visible channel. This is left as zero for the IR channel.
    nphangle = len(exop.phangle)
    pha = np.zeros((nprofiles,nchan,nlevels-1,nphangle), dtype=np.float64)
    pha[0,0,74,:] = exop.phasefn[0,:]     # phase fn: profile 1, channel 1, layer 75
    pha[1,0,74,:] = exop.phasefn[0,:]     # phase fn: profile 2, channel 1, layer 75

    # Calculate Legendre coefficients for visible channel as DOM is being used. We only need
    # DomNstreams coefficients (note that this requires a DomNstreams+1 sized array: the zeroth
    # coefficient is always 1.) These are left as zero for the IR channel.
    nmom = seviriRttov.Options.DomNstreams
    legcoef = np.zeros((nprofiles,nchan,nlevels-1,nmom+1), dtype=np.float64)
    lc = seviriRttov.calcLegcoef(exop.phangle, exop.phasefn[0,:], nmom)
    legcoef[0,0,74,:] = lc[:]             # Leg. coefs: profile 1, channel 1, layer 75
    legcoef[1,0,74,:] = lc[:]             # Leg. coefs: profile 2, channel 1, layer 75

    # Now assign the optical property arrays to the Rttov object
    seviriRttov.CldAsb     = asb
    seviriRttov.CldPhangle = exop.phangle
    seviriRttov.CldPha     = pha
    seviriRttov.CldLegcoef = legcoef


    # Associate the profiles with the Rttov instance
    seviriRttov.Profiles = myProfiles

    # ------------------------------------------------------------------------
    # Call RTTOV
    # ------------------------------------------------------------------------

    # In this example we let the Rttov class automatically set calcemis/calcrefl
    # to true for simplicity. You can supply emissivity/BRDF values or use the
    # atlases in exactly the same way as for other simulation types.

    # Call the RTTOV direct model: note the channel list must be consistent with the
    # channels for which optical properties are specified.
    try:
        seviriRttov.runDirect(chan_list)
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)

    # ------------------------------------------------------------------------
    # Print out some of the output
    # ------------------------------------------------------------------------

    for p in range(nprofiles):
        print("Profile {:d}:".format(p))
        print("  Visible channel cloudy and clear reflectances:")
        print("  {:f}  {:f}".format(seviriRttov.Refl[p,0], seviriRttov.ReflClear[p,0]))
        print("  IR channel cloudy and clear BTs:")
        print("  {:f}  {:f}".format(seviriRttov.Bt[p,1], seviriRttov.BtClear[p,1]))
        print

    # ------------------------------------------------------------------------
    # Deallocate memory
    # ------------------------------------------------------------------------

    # Because of Python's garbage collector, there should be no need to
    # explicitly deallocate memory
