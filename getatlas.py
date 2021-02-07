#coding:utf-8

import pyrttov

def GetBRDF(AtlasPath, month, myrttov):
    brdfAtlas = pyrttov.Atlas()
    brdfAtlas.AtlasPath = AtlasPath
    brdfAtlas.loadBrdfAtlas(month, myrttov) # Supply Rttov object to enable single-instrument initialisation
    brdfAtlas.IncSea = False                # Do not use BRDF atlas for sea surface types

    return brdfAtlas

def GetIREmis(AtlasPath, month):

    irAtlas = pyrttov.Atlas()
    irAtlas.AtlasPath = AtlasPath
    irAtlas.loadIrEmisAtlas(month, ang_corr=True) # Include angular correction, but do not initialise for single-instrument

    return irAtlas

def GetMWEmis(AtlasPath, month):
    # TELSEM2 atlas does not require an Rttov object to initialise
    mwAtlas = pyrttov.Atlas()
    mwAtlas.AtlasPath = AtlasPath
    mwAtlas.loadMwEmisAtlas(month)

    return mwAtlas



