#coding:utf-8
import os
import xarray as xr
import numpy as np
from config import *

from scipy import interpolate

from ncpro import lb_writenc,lb_readnc
from rttov_ir_clr_fwd import SimulateIR




class CSimulateT799(object):

    def __init__(self, path, nowdate, height, width,
                 satzen, satazi, sunzen, sunazi,
                 lat, lon, elev, surftype, year, month, day):

        self.Height = height
        self.Width = width

        # 将经度范围调整为0 ~ 360
        lon[(lon>=-180.0) & (lon<0)] += 360

        self.pathin = path
        self.nowdate = nowdate

        self.satz = satzen
        self.sata = satazi
        self.sunz = sunzen
        self.suna = sunazi
        self.lat = lat
        self.lon = lon
        self.dem = elev * 0.001     # m-->Km
        self.year = year
        self.month = month
        self.day = day

        '''
        0 = Shallow Ocean (Ocean < 5km from coast or < 50m deep).
        1 = Lamd (not anything else).
        2 = Ocean Coastlines and Lake Shorelines.
        3 = Shallow Inland Water (Inland Water < 5km from shore or < 50m deep).
        4 = Ephemeral (intermittent) Water.
        5 = Deep Inland Water ( Inland water > 5km from shoreline and > 50m deep).
        6 = Moderate or Continental Ocean (Ocean > 5km from coast and > 50m deep and <500m deep).
        7 = Deep Ocean (Ocean > 500m deep). 
        '''

        self.surftype = np.full((2, self.Height, self.Width), fill_value=0, dtype=np.int32)
        self.surftype[0, surftype!=1] = 1  # Surface type (land = 0, sea = 1, seaice = 2)
        self.surftype[1, surftype!=1] = 1  # Water type (fresh = 0, ocean = 1)


        # 对输入有效性进行检查
        self.validflag = self.DataValidCheck()
        self.nprofs = np.sum(self.validflag)

        # 读取NWP数据
        self.GetProfData(path, nowdate)

    def rttov(self, nchans, channelID, FileCoef,
              emis=None, brdf=None,
              emispath=None, brdfpath=None):

        # 视场匹配
        row, col = self.ViewMatch()
        self.nchans = nchans

        if emis is not None :
            inemis = np.transpose(emis[:, row, col])
        else:
            inemis = None

        if brdf is not None :
            inbrdf= np.full((len(row), nchans), fill_value=-1.0, dtype=np.float32)
            inbrdf[:, :6] = np.transpose(brdf[:, row, col])
        else:
            inbrdf = None

        self.myrttov = SimulateIR(
            nchans, channelID,
            self.nprofs, self.NWP_LEVEL,
            np.transpose(self.nwp_pre[:, row, col]),
            np.transpose(self.nwp_t[:, row, col]),
            np.transpose(self.nwp_q[:, row, col]),
            self.nwp_t2m[row, col], self.nwp_q2m[row, col],
            self.nwp_ts[row, col], self.nwp_sp[row, col],
            self.nwp_u10m[row, col], self.nwp_v10m[row, col],
            self.satz[self.validflag], self.sata[self.validflag],
            self.sunz[self.validflag], self.suna[self.validflag],
            self.lat[self.validflag], self.lon[self.validflag],
            self.dem[self.validflag], self.surftype[:, self.validflag],
            self.year, self.month,
            coeffilename=FileCoef,
            inemis=inemis, inbrdf=inbrdf,
            emispath=emispath, brdfpath=brdfpath)

        return self.myrttov

        # myprof = Initprofiles(self.nprofs, self.NWP_LEVEL,
        #                       np.transpose(self.nwp_pre[:, row, col]),
        #                       np.transpose(self.nwp_t[:, row, col]),
        #                       np.transpose(self.nwp_q[:, row, col]),
        #                       self.nwp_t2m[row, col],
        #                       self.nwp_q2m[row, col],
        #                       self.nwp_ts[row, col],
        #                       self.nwp_sp[row, col],
        #                       self.nwp_u10m[row, col],
        #                       self.nwp_v10m[row, col],
        #                       # self.fy4a_nom_p[self.validflag, :],
        #                       # self.fy4a_nom_t[self.validflag, :],
        #                       # self.fy4a_nom_q[self.validflag, :],
        #                       # self.fy4a_nom_t2m[self.validflag],
        #                       # self.fy4a_nom_q2m[self.validflag],
        #                       # self.fy4a_nom_ts[self.validflag],
        #                       # self.fy4a_nom_sp[self.validflag],
        #                       # self.fy4a_nom_u10m[self.validflag],
        #                       # self.fy4a_nom_v10m[self.validflag],
        #                       self.satz[self.validflag],
        #                       self.sata[self.validflag],
        #                       self.sunz[self.validflag],
        #                       self.suna[self.validflag],
        #                       self.lat[self.validflag],
        #                       self.lon[self.validflag],
        #                       self.dem[self.validflag],
        #                       self.surftype[:, self.validflag],
        #                       self.year,
        #                       self.month
        #                       )


        # self.myrttov = SimulateIR_clr(nchans, channelID, self.nprofs,
        #                self.NWP_LEVEL, myprof, self.month, FileCoef,
        #                inemis=inemis, inbrdf=inbrdf,
        #                emispath=emispath, brdfpath=brdfpath)


    def GetBT(self):
        # np.save('btref.npy', self.myrttov.BtRefl)
        bt = np.full((self.nchans, self.Height, self.Width), dtype=np.float32, fill_value=-999.0)
        for ichan in range(self.nchans) :
            bt[ichan, self.validflag] = self.myrttov.BtRefl[:, ichan]

        # rbt = np.load('btref.npy')
        # bt = np.full((self.nchans, self.Height, self.Width), dtype=np.float32, fill_value=-999.0)
        # for ichan in range(self.nchans) :
        #     bt[ichan, self.validflag] = rbt[:, ichan]


        return bt

    def Interp2View(self, y, x, z, y1, x1):
        '''
        插值到视场范围
        :param x:
        :param y:
        :param z:
        :param x1:
        :param y1:
        :return:
        '''
        # print(len(x))
        # print(len(y))
        # print(z.shape)

        fit = interpolate.interp2d(x, y, z)

        return fit(x1, y1)

    def ViewMatch(self):

        row = np.array((90.0 - self.lat[self.validflag]) / self.nwp_resolution, dtype=np.int32)
        col = np.array((self.lon[self.validflag] - 0.) / self.nwp_resolution, dtype=np.int32)

        # 此模块暂不放开，测试程序错误情况时的经纬度 #?????????
        # row[(row < 0) ] = 0
        # col[(col < 0) ] = 0
        #
        # row[(row >= self.NWP_LINE)] = self.NWP_LINE - 1
        # col[(col >= self.NWP_PIXEL)] = self.NWP_PIXEL - 1
        return row, col

        # fy4a_lat = np.arange(81.222336,  -81.222336, -0.08)
        # fy4a_lon = np.arange(23.422451, -174.02245 + 360, 0.08)
        #
        # LINE_4000 = lb_readnc('./parm/LATLON_TO_PIXLINE_4000M.HDF', 'LINE_4000')
        # PIX_4000 = lb_readnc('./parm/LATLON_TO_PIXLINE_4000M.HDF', 'PIX_4000')
        # flag = (LINE_4000 !=  65535.0) & (PIX_4000 != 65535.0)
        # row = LINE_4000[flag]
        # col = PIX_4000[flag]
        #
        #
        # self.fy4a_nom_p    = np.full((self.Height, self.Width, self.NWP_LEVEL), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_t    = np.full((self.Height, self.Width, self.NWP_LEVEL), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_q    = np.full((self.Height, self.Width, self.NWP_LEVEL), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_t2m  = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_q2m  = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_ts   = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_sp   = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_u10m = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        # self.fy4a_nom_v10m = np.full((self.Height, self.Width), fill_value=-999.0, dtype=np.float64)
        #
        # for ilev in range(self.NWP_LEVEL) :
        #     self.fy4a_nom_p[row, col, ilev] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_pre[ilev,:,:], fy4a_lat, fy4a_lon)[flag]
        #     self.fy4a_nom_t[row, col, ilev] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_q[ilev,:,:], fy4a_lat, fy4a_lon)[flag]
        #     self.fy4a_nom_q[row, col, ilev] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_q[ilev,:,:], fy4a_lat, fy4a_lon)[flag]
        #
        # self.fy4a_nom_t2m [row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_t2m,  fy4a_lat, fy4a_lon)[flag]
        # self.fy4a_nom_q2m [row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_q2m,  fy4a_lat, fy4a_lon)[flag]
        # self.fy4a_nom_ts  [row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_ts,   fy4a_lat, fy4a_lon)[flag]
        # self.fy4a_nom_sp  [row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_sp,   fy4a_lat, fy4a_lon)[flag]
        # self.fy4a_nom_u10m[row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_u10m, fy4a_lat, fy4a_lon)[flag]
        # self.fy4a_nom_v10m[row, col] = self.Interp2View(self.nwp_lat, self.nwp_lon, self.nwp_v10m, fy4a_lat, fy4a_lon)[flag]

        # lb_writenc('nwp_fy4a_nom.nc', 'lat', range(2748), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'lon', range(2748), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 't', self.fy4a_nom_t,dimension=('lat', 'lon', 'lev'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'q', self.fy4a_nom_q,dimension=('lat', 'lon', 'lev'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 't2m', self.fy4a_nom_t2m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'q2m', self.fy4a_nom_q2m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'ts', self.fy4a_nom_ts,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'sp', self.fy4a_nom_sp,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'u10m', self.fy4a_nom_u10m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp_fy4a_nom.nc', 'v10m', self.fy4a_nom_v10m,dimension=('lat', 'lon'), overwrite=0)
        # exit()

    def DataValidCheck(self):

        return (self.sata >= 0) & (self.sata <= 360.0) \
               & (self.satz >= 0) & (self.satz <= 75.0) \
               & (self.sunz >= 0) & (self.sunz <= 180.0) \
               & (self.suna >= 0) & (self.suna <= 360.0) \
               & (self.lat >= -90.0) & (self.lat <= 90.0) \
               & (self.lon >= 0.0) & (self.lon <= 360.0)

    def GetProfData(self, path, nowdate):

        self.nwp_sp = self.readT799(path, 'KTPSG', nowdate) * 0.01 # Pa --> hPa
        nwp_rh2m = self.readT799(path, 'KTQ2G', nowdate)
        nwp_rh = self.readT799(path, 'KTRHG', nowdate)[::-1,:,:]

        nwp_rh2m[nwp_rh2m>100] = 100
        nwp_rh2m[nwp_rh2m<0.0] = 0.01

        nwp_rh[nwp_rh>100] = 100
        nwp_rh[nwp_rh<0.0] = 0.01

        self.nwp_t2m = self.readT799(path, 'KTT2G', nowdate)
        self.nwp_ts = self.readT799(path, 'KTTSG', nowdate)
        self.nwp_t = self.readT799(path, 'KTTTG', nowdate)[::-1,:,:]
        self.nwp_u10m = self.readT799(path, 'KTUTG', nowdate)
        self.nwp_v10m = self.readT799(path, 'KTVTG', nowdate)

        self.nwp_lat = self.readT799(path, 'KTTTG', nowdate, 'latitude')
        self.nwp_lon = self.readT799(path, 'KTTTG', nowdate, 'longitude')
        pre = self.readT799(path, 'KTTTG', nowdate, 'isobaricInhPa')[::-1]

        self.NWP_LINE = len(self.nwp_lat)
        self.NWP_PIXEL = len(self.nwp_lon)
        self.NWP_LEVEL = len(pre)
        self.nwp_resolution = 360.0/ self.NWP_PIXEL

        self.nwp_pre = np.full_like(nwp_rh, fill_value=-999.)
        for ilev in range(self.NWP_LEVEL) :
            self.nwp_pre[ilev, :, :] = pre[ilev]

        # 根据温度和相对湿度，计算比湿 （kg/kg）
        self.nwp_q = self.ESA(self.nwp_t, nwp_rh, self.nwp_pre)
        self.nwp_q2m = self.ESA(self.nwp_t2m, nwp_rh2m, self.nwp_sp)

        # self.nwp_q[self.nwp_q<1.0E-6] = 1.0E-6
        # self.nwp_q2m[self.nwp_q2m<1.0E-6] = 1.0E-6


        self.nwp_q *= 1e6*28.966/18.0 # kg/kg ==> ppmv
        self.nwp_q2m *= 1e6*28.966/18.0 # kg/kg ==> ppmv


        # lb_writenc('nwp.nc', 'lev', pre, overwrite=1)
        # lb_writenc('nwp.nc', 'lat', self.nwp_lat, overwrite=0)
        # lb_writenc('nwp.nc', 'lon', self.nwp_lon, overwrite=0)
        # lb_writenc('nwp.nc', 't', self.nwp_t,dimension=('lev', 'lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'q', self.nwp_q,dimension=('lev', 'lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'q2m', self.nwp_q2m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 't2m', self.nwp_t2m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'ts', self.nwp_ts,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'sp', self.nwp_sp,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'u10m', self.nwp_u10m,dimension=('lat', 'lon'), overwrite=0)
        # lb_writenc('nwp.nc', 'v10m', self.nwp_v10m,dimension=('lat', 'lon'), overwrite=0)
        # exit()

    def readT799(self, pathin, key, nowdate, sdsname=None):

        if key in T799 :
            filename = os.path.join(pathin, '{key}{pretime}'.format(
                key=key,
                pretime=nowdate
            ))
            if os.path.isfile(filename) :
                ds = xr.open_dataset(filename, engine="cfgrib")
                # ds_grib = xr.open_dataarray(filename, engine="cfgrib")
                # ds.to_netcdf(filename+'.nc')
                if sdsname is None :
                    for item in ds.data_vars.keys():
                        a = ds.get(item)
                        print('read %s  %s success...' %(filename, item))
                    return a.data
                else:
                    a = ds.get(sdsname)
                    print('read %s  %s success...' %(filename, sdsname))
                    return a.data
            else:
                print('%s is not exist!!!' %(filename))
                return None


    def ESA(self, temp, rh, presure):
        '''
            计算比湿
        :return: kg/kg
        '''
        flag = temp < 273.16
        q = temp.copy()

        q[flag] = 6.11 * np.exp(17.67 * (temp[flag] - 273.16) / (temp[flag] - 29.66))
        q[~flag] = 6.11 * np.exp(17.269 * (temp[~flag] - 273.16) / (temp[~flag] - 35.86))

        fEp = q * rh * 0.01

        return  (fEp / (presure - 0.378 * fEp)) * 0.622

def run(nwptime, geoname, outname, albedoname=None):
    import time

    t1 = time.time()
    if albedoname is not None :
        brdf = lb_readnc(albedoname, 'brdf')


    satzen = lb_readnc(geoname, 'NOMSatelliteZenith')
    satazi = lb_readnc(geoname, 'NOMSatelliteAzimuth')
    sunzen = lb_readnc(geoname, 'NOMSunZenith')
    sunazi = lb_readnc(geoname, 'NOMSunAzimuth')
    lat = lb_readnc('/product/mnt/test/parm/FY4A_OBI_4000M_NOM_LATLON.HDF', 'Lat')
    lon = lb_readnc('/product/mnt/test/parm/FY4A_OBI_4000M_NOM_LATLON.HDF', 'Lon')
    elev = lb_readnc('/product/mnt/test/parm/IFL_FY4A_AGRIX_DEM_4000M.HDF', 'DEM')
    LandSeaMask = lb_readnc('/product/mnt/test/parm/IFL_FY4A_AGRIX_LMK_4000M.HDF', 'LandSeaMask')

    year = 2020
    month = 8
    day = 17
    msimul = CSimulateT799(path,nwptime, 2748, 2748,
                           satzen, satazi, sunzen, sunazi,
                           lat, lon, elev, LandSeaMask, year, month, day)

    rttov_installdir = r'/product/work/reflib/rttov13/'
    emispath = '{}/{}/{}'.format('/product/work/parm', "emis_data",'IR')
    brdfpath = '{}/{}'.format('/product/work/parm', "brdf_data")
    FileCoef = '{}/{}'.format(rttov_installdir,
                              "rtcoef_rttov13/rttov13pred54L/rtcoef_fy4_1_agri_o3co2.dat")
    nchans = 13
    channels = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

    msimul.rttov(nchans, channels, FileCoef,
                 brdf=brdf, emispath=emispath, brdfpath=brdfpath)

    bt = msimul.GetBT()

    lb_writenc(outname, 'bt', bt, dimension=('z', 'y', 'x'), overwrite=1)

    t2 = time.time()
    print('cost %.2fs' %(t2 - t1))
    exit()

if __name__ == '__main__':
    path = '/product/mnt/test/data/20200817'

    import glob
    import datetime
    filelist = glob.glob(os.path.join(r'/product/mnt/data/simulate/GEO', '*L1*GEO*.HDF'))
    filelist.sort()
    for geoname in filelist :
        basename = os.path.basename(geoname)
        namelist = basename.split('_')
        nowtime = datetime.datetime.strptime(namelist[9], '%Y%m%d%H%M%S')
        albedoname = r'/product/mnt/data/albedo/fy4a_nom_%s.nc' %(namelist[9])

        hours = nowtime.hour
        nwptime ='%s%02d.%03d' %('20200817' , int(hours/12) * 12, int(np.mod(hours, 12) / 3) * 3)
        outname = os.path.join(r'/product/mnt/data/simulate', 'fy4a_nom_brdf_%s.nc' %(nowtime.strftime('%Y%m%d%H%M%S')))
        run(nwptime, geoname, outname, albedoname)