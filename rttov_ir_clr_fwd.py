#coding:utf-8
import sys
sys.path.append(r'/product/work/reflib/rttov13/wrapper')
import pyrttov
import numpy as np
import time
from getatlas import GetIREmis, GetBRDF

FILLVALUE = 0.0


'''
In particular the input
cloud and hydrometeor arrays are defined on nlevels (unlike the case for visible/IR scattering where
they are on nlayers) and the pressure half-levels profile has size (nlevels+1). 
'''
def Initprofiles(nprofiles, nlevels,
                 pressure, t, q,
                 t2m, q2m,
                 ts, ps,
                 u10m, v10m,
                 satz, sata,
                 sunz, suna,
                 lat, lon,
                 dem, landcover,
                 year, month,
                 ctp=None, cfrac=None,
                 clw_scheme=None, clwde_param=None,
                 Be=None, cosbk=None,
                 isclear = True,
                 O3=None, CO2=None):
    # 初始化廓线
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    # 廓线对象包含气压、温度、湿度、CO2
    myProfiles.GasUnits = 2  # ppmv over moist air  ppmv
    myProfiles.P = pressure  # [nprof, nlevels]     hPa
    myProfiles.T = t  # [nprof, nlevels]            K
    myProfiles.Q = q  # [nprof, nlevels]            ppmv
    if CO2 is not None :
        myProfiles.CO2 = CO2 # [nprof, nlevels]     ppmv
    if O3 is not  None :
        myProfiles.O3 = O3   # [nprof, nlevels]     ppmv
    #--------------------------------------------------------------------------------
    angles = np.full((nprofiles, 4), dtype=np.float64, fill_value=FILLVALUE)
    angles[:, 0] = satz
    angles[:, 1] = sata
    angles[:, 2] = sunz
    angles[:, 3] = suna
    myProfiles.Angles = angles # angles[4][nprofiles]: satzen, satazi, sunzen, sunazi
    #--------------------------------------------------------------------------------
    s2m = np.full((nprofiles, 6), dtype=np.float64, fill_value=FILLVALUE)
    s2m[:, 0] = ps
    s2m[:, 1] = t2m
    s2m[:, 2] = q2m
    s2m[:, 3] = u10m
    s2m[:, 4] = v10m
    s2m[:, 5] = 100000.
    myProfiles.S2m = s2m       # 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
    #--------------------------------------------------------------------------------
    skin = np.full((nprofiles, 9), dtype=np.float64, fill_value=0)
    skin[:, 0] = ts
    skin[:, 4] = 3
    skin[:, 5] = 5
    skin[:, 6] = 15
    skin[:, 7] = 0.1
    skin[:, 8] = 0.3

    myProfiles.Skin = skin       # Skin T, salinity, snow_fraction, foam_fraction, fastem_coefs(1:5).
    #--------------------------------------------------------------------------------
    surftype = np.full((nprofiles, 2), dtype=np.int32, fill_value=0)
    # surftype[:, 0] = 1          #?????????????
    # surftype[:, 1] = 0          #?????????????
    surftype[:, 0] = landcover[0,:]
    surftype[:, 1] = landcover[1,:]
    myProfiles.SurfType = surftype # surftype [2][nprofiles]: surftype, watertype
    #--------------------------------------------------------------------------------
    surfgeom = np.full((nprofiles, 3), dtype=np.float64, fill_value=FILLVALUE)
    surfgeom[:, 0] = lat
    surfgeom[:, 1] = lon
    surfgeom[:, 2] = dem
    myProfiles.SurfGeom = surfgeom # surfgeom [3][nprofiles]: lat, lon, elev
    #--------------------------------------------------------------------------------
    datetimes = np.full((nprofiles, 6), dtype=np.int32, fill_value=0)
    datetimes[:, 0] = year
    datetimes[:, 1] = month
    datetimes[:, 2] = 1
    datetimes[:, 3] = 0
    datetimes[:, 4] = 0
    datetimes[:, 5] = 0
    myProfiles.DateTimes = datetimes # Year, month, day, hour, minute, second

    # Optional
    # myProfiles.SimpleCloud  # real  [nprofiles][2]  (ctp, cfraction) per profile
    # myProfiles.ClwScheme   #  Integer [nprofiles][2]  Visible/IR (clw_scheme, clwde_param) per profile
    # myProfiles.IceCloud    #  Integer [nprofiles][2]  (ice_scheme, icede_param) per profile
    # myProfiles.Zeeman      #  Real [nprofiles][2]  (Be, cosbk) per profile

    # myProfiles.Cfrac    #cfrac  (cloud fraction)
    # myProfiles.Cirr =   # ciw  (cloud ice water)
    # myProfiles.Inso  = # aer_inso (insoluble aerosol)

    return myProfiles


def SimulateIR_clr(nchans, channels, nprofiles, nlevels,
                   myProfiles,
                   month, coeffilename=None,
                   inemis=None, inbrdf=None,
                   emispath=None, brdfpath=None,
                   Nthreads= 10,
                   btFlag=False,
                   radianceFlag=False,
                   tauTotalFlag=False):

    # 初始化RTTOV模型
    myrttov = pyrttov.Rttov()

    # 初始化系数文件
    myrttov.FileCoef = coeffilename
    # myrttov.FileScaer
    # 初始化模式设置
    # Use the RTTOV interpolator
    myrttov.Options.AddInterp = True
    myrttov.Options.InterpMode = 2
    # Turn on solar radiation
    myrttov.Options.AddSolar = True
    # Turn on verbose wrapper output
    myrttov.Options.VerboseWrapper = True
    # Take advantage of multiple threads if RTTOV was compiled with OpenMP
    myrttov.Options.Nthreads = Nthreads

    myrttov.Options.IrSeaEmisModel = 1
    myrttov.Options.IrScattModel = 1
    myrttov.Options.VisScattModel = 1

    myrttov.Options.DoCheckinput = True
    myrttov.Options.Verbose = False

    # 加载通道信息
    try:
        myrttov.loadInst(channels)
    except pyrttov.RttovError as e:
        sys.stderr.write("Error loading instrument(s): {!s}".format(e))
        sys.exit(1)

    # 初始化廓线数据
    myrttov.Profiles = myProfiles

    # ------------------------------------------------------------------------
    # Call RTTOV
    # ------------------------------------------------------------------------

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
    # calcemis/calcrefl TRUE - see RTTOV user guide)

    # 初始化地表发射率，其维度为array（4, nprofiles, nchans）,
    # emissivity (index 0), BRDF (index 1), diffuse reflectance (index 2), and
    # specularity (index 3) for all channels and profiles being simulated
    surfemisrefl_ir = np.full((4, nprofiles, nchans), dtype=np.float64, fill_value=-1)
    myrttov.SurfEmisRefl = surfemisrefl_ir

    seaflag = myProfiles.SurfType[:,0] == 1

    if inemis is None:
        try:
            # 计算发射率
            irEmisAtlas = GetIREmis(emispath, month)

            # Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
            # *loaded* channels which is what is required
            surfemisrefl_ir[0,:,:] = irEmisAtlas.getEmisBrdf(myrttov)
        except pyrttov.RttovError as e:
            # If there was an error the emissivities/BRDFs will not have been modified so it
            # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
            sys.stderr.write("Error calling atlas: {!s}".format(e))
    else:
        surfemisrefl_ir[0,:,:] = inemis

    if inbrdf is None :
        try:
            # 计算BRDF
            brdfAtlas = GetBRDF(brdfpath, month, myrttov)

            surfemisrefl_ir[1,:,:] = brdfAtlas.getEmisBrdf(myrttov)
        except pyrttov.RttovError as e:
            # If there was an error the emissivities/BRDFs will not have been modified so it
            # is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
            sys.stderr.write("Error calling atlas: {!s}".format(e))
    else:
        surfemisrefl_ir[1,:,:] = inbrdf

    surfemisrefl_ir[:, seaflag, :] = -1

    # Call the RTTOV direct model for each instrument:
    # no arguments are supplied to runDirect so all loaded channels are
    # simulated
    try:
        t1 = time.time()
        myrttov.runDirect(channels)
        t2 = time.time()
        print('runDirect cost %.2fs' %(t2 - t1))
    except pyrttov.RttovError as e:
        sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
        sys.exit(1)

    # Call the RTTOV K model for each instrument:
    # Jacobian
    # try:
    #     myrttov.runK(channels)
    # except pyrttov.RttovError as e:
    #     sys.stderr.write("Error running RTTOV direct model: {!s}".format(e))
    #     sys.exit(1)

    # bt = myrttov.BtRefl
    # radiance = myrttov.Rads
    # tau = myrttov.TauTotal

    # print(bt)
    # print('*'*100)
    # print(myrttov.RadTotal)
    #
    # print('*'*100)

    return myrttov


