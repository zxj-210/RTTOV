#coding:utf-8
'''
涉及可见光、红外、微波、雷达
区分晴空、云区

'''

import sys
import os
sys.path.append(r'/home/htht/rttov/rttov13/wrapper')
import pyrttov
# os.environ[''] = r'/home/htht/rttov/rttov13/wrapper'
from example import example_data as ex
import numpy as np

from rttov_ir_clr_fwd import SimulateIR_clr, Initprofiles

# from GetDataT799 import GetDataFromT799

def expand2nprofiles(n, nprof):
    '''
    获取测试数据
    :param n:
    :param nprof:
    :return:
    '''
    """Transform 1D array to a [nprof, nlevels] array"""
    outp = np.empty((nprof, len(n)), dtype=n.dtype)
    for i in range(nprof):
        outp[i, :] = n[:]
    return outp


def test():

    rttov_installdir = r'/home/htht/rttov/rttov13/'
    # Declare an instance of Profiles
    nlevels = len(ex.p_ex)

    nprofiles = 2
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # Associate the profiles and other data from example_data.h with myProfiles
    # Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and
    # are omitted here


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
    return myProfiles
    # ------------------------------------------------------------------------
    # Set up Rttov instances for each instrument
    # ------------------------------------------------------------------------

    # Create three Rttov objects for three instruments
    seviriRttov = pyrttov.Rttov()

    # For HIRS and MHS we will read all channels, but we will read a subset
    # for SEVIRI
    nchan_seviri = 10


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



    # Load the instruments: for HIRS and MHS do not supply a channel list and
    # so read all channels
    try:
        seviriRttov.loadInst(chan_list_seviri)

    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    # Associate the profiles with each Rttov instance
    seviriRttov.Profiles = myProfiles


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



    # Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    surfemisrefl_seviri = np.zeros((4,nprofiles,nchan_seviri), dtype=np.float64)

    seviriRttov.SurfEmisRefl = surfemisrefl_seviri
  #-----------------------------------------------

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
    # calcemis/calcrefl TRUE - see RTTOV user guide)

    surfemisrefl_seviri[:,:,:] = -1.

    # Call emissivity and BRDF atlases
    try:
        # Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
        # *loaded* channels which is what is required
        surfemisrefl_seviri[0,:,:] = irAtlas.getEmisBrdf(seviriRttov)
        surfemisrefl_seviri[1,:,:] = brdfAtlas.getEmisBrdf(seviriRttov)

    except pyrttov.RttovError as e:
        # If there was an error the emissivities/BRDFs will not have been modified so it
        # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        sys.stderr.write("Error calling atlas: {!s}".format(e))

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are
    # simulated
    try:
        seviriRttov.runDirect()

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



    # Because of Python's garbage collector, there should be no need to
    # explicitly deallocate memory


if __name__ == '__main__':

    # proffilename = r'KTGHG2019080100.000'
    # GetDataFromT799(proffilename)
    # exit()

    rttov_installdir = r'/home/htht/rttov/rttov13/'
    emispath = '{}/{}/{}'.format(rttov_installdir, "emis_data",'IR')
    brdfpath = '{}/{}'.format(rttov_installdir, "brdf_data")
    FileCoef = '{}/{}'.format(rttov_installdir,
                   "rtcoef_rttov13/rttov13pred54L/rtcoef_fy4_1_agri_o3co2.dat")
    nchans = 13
    channels = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

    # 廓线条数
    nprofiles = 2
    # 廓线层数
    nlevels = 101

    P = expand2nprofiles(ex.p_ex, nprofiles)
    T = expand2nprofiles(ex.t_ex, nprofiles)
    # Modify the temperature of the second profile
    T[1, :] += 2
    Q = expand2nprofiles(ex.q_ex, nprofiles)
    CO2 = expand2nprofiles(ex.co2_ex, nprofiles)

    satzen  = ex.angles[:,0]
    satazi  = ex.angles[:,1]
    sunzen  = ex.angles[:,2]
    sunazi  = ex.angles[:,3]

    ps = ex.s2m[:, 0]
    t2m = ex.s2m[:, 1]
    q2m = ex.s2m[:, 2]
    u10m = ex.s2m[:, 3]
    v10m = ex.s2m[:, 4]

    ts = ex.skin[:, 0]

    surftype = ex.surftype[:, 0]
    watertype = ex.surftype[:,1]

    lat = ex.surfgeom[:, 0]
    lon = ex.surfgeom[:, 1]
    elev = ex.surfgeom[:, 2]

    years = ex.datetimes[:, 0]
    months = ex.datetimes[:, 1]
    days = ex.datetimes[:, 2]
    hours = ex.datetimes[:,3]
    minuts = ex.datetimes[:, 4]
    seconds = ex.datetimes[:,5]

    myprof = Initprofiles(nprofiles, nlevels, P, T, Q, t2m, q2m, ts, ps, u10m, v10m,
                          satzen, satazi, sunzen, sunazi, lat, lon, elev, surftype, 2017, month, CO2=CO2)
    # myprof = test()

    bt = SimulateIR_clr(nchans, channels, nprofiles, nlevels, myprof, 8,
                   FileCoef, emispath=emispath, brdfpath=brdfpath, btFlag=True)

    print(bt)

