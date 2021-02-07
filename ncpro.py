#coding:utf-8

import os
import netCDF4
import numpy as np


def lb_readnc(filename, sdsname,  dictsdsinfo=None):

    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return None
    else:
        with netCDF4.Dataset(filename, 'r') as fp :
            if not sdsname in fp.variables:
                print('%s is not a dataset in %s' % (sdsname, filename))
                return None
            else:
                dsetid = fp.variables[sdsname]
                data = fp.variables[sdsname][:]

                if not dictsdsinfo is None :
                    for key in dsetid.ncattrs() :
                        dictsdsinfo.update({key: dsetid.getncattr(key)})

                    return data, dictsdsinfo

        return data

def lb_readnc_fileinfo(filename) :

    dictfileinfo = {}
    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return dictfileinfo
    else:
        with netCDF4.Dataset(filename, 'r') as fp :

            for key in fp.ncattrs():
                dictfileinfo.update({key: fp.getncattr(key)})

            return  dictfileinfo

def lb_readnc_sdsinfo(filename, sdsname) :

    dictsdsinfo = {}
    if not os.path.isfile(filename) :
        print('%s is not exist, please check it...' %(filename))
        return dictsdsinfo
    else:
        with netCDF4.Dataset(filename, 'r') as fp :

            if not sdsname in fp.variables:
                print('%s is not a dataset in %s' %(sdsname, filename))
                return dictsdsinfo
            else:
                dsetid = fp.variables[sdsname]
                for key in dsetid.ncattrs():
                    dictsdsinfo.update({key: dsetid.getncattr(key)})

                return dictsdsinfo


def lb_writenc(filename, sdsname, srcdata, dimension=None, overwrite=1, complevel=9,  dictsdsinfo=None):
    '''
    写NC文件
    :param filename: 输出文件名
    :param sdsname: 数据集名
    :param srcdata: 数据
    :param dimension: tuple, 关联维度名，
    :param overwrite:
    :param complevel:
    :param dictfileinfo:
    :param dictsdsinfo:
    :return:
    '''
    data = np.array(srcdata)

    if overwrite :
        fp = netCDF4.Dataset(filename, 'w')
    else:
        fp = netCDF4.Dataset(filename, 'r+')

    dsetid = None
    # 如果没有对应的关联维度，则创建该维度
    if dimension is None :
        fp.createDimension(sdsname, data.shape[0])
        dsetid = fp.createVariable(sdsname, data.dtype, dimensions=sdsname,
                                  zlib=True, complevel=9)
        fp.variables[sdsname][:] = data
    else:
        # 判断数据集维度是否存在，如果不存在，则创建该维度，并赋值为0
        k = 0
        for item in dimension :
            if not item in fp.dimensions :
                fp.createDimension(item, data.shape[k])
                dsetid = fp.createVariable(item, data.dtype, dimensions=item, zlib=True, complevel=9)
                fp.variables[item][:] = range(data.shape[k])
            k = k + 1
        dsetid = fp.createVariable(sdsname, data.dtype, dimensions=dimension, zlib=True, complevel=complevel)
        fp.variables[sdsname][:] = data

    if (not dictsdsinfo is None) and (not dsetid is None) :
        for key in dictsdsinfo :
            dsetid.setncattr(key, dictsdsinfo[key])

    fp.close()


def lb_writenc_fileinfo(filename, dictfileinfo,overwrite=1):

    try:
        if overwrite :
            fp = netCDF4.Dataset(filename, 'w')
        else:
            if os.path.isfile(filename) :
                fp = netCDF4.Dataset(filename, 'r+')
            else:
                print('%s is not exist, will create it...' %(filename))
                fp = netCDF4.Dataset(filename, 'w')

        if not dictfileinfo is None :
            for key in dictfileinfo:
                fp.setncattr(key, dictfileinfo[key])
        fp.close()
        return True

    except BaseException as e :
        print('write %s file information error!!!' %(filename))
        return False


def lb_writenc_sdsinfo(filename, sdsname,  dictsdsinfo, overwrite=1):
    try:
        if overwrite:
            fp = netCDF4.Dataset(filename, 'w')
        else:
            if os.path.isfile(filename):
                fp = netCDF4.Dataset(filename, 'r+')
            else:
                print('%s is not exist, will create it...' % (filename))
                fp = netCDF4.Dataset(filename, 'w')
        if not sdsname in fp.variables:
            print('%s is not a dataset in %s' %(sdsname, filename))
            fp.close()
            return False
        else:
            dsetid = fp.variables[sdsname]
            for key in dictsdsinfo:
                dsetid.setncattr(key, dictsdsinfo[key])
        fp.close()
        return True
    except BaseException as e:
        print('write %s dataset %s information error!!!' % (filename, sdsname))
        return False



# with netCDF4.Dataset("in.nc") as src, netCDF4.Dataset("out.nc", "w") as dst:
#     # copy global attributes all at once via dictionary
#     dst.setncatts(src.__dict__)
#     # copy dimensions
#     for name, dimension in src.dimensions.items():
#         dst.createDimension(
#             name, (len(dimension) if not dimension.isunlimited() else None))
#     # copy all file data except for the excluded
#     for name, variable in src.variables.items():
#         if name not in toexclude:
#             x = dst.createVariable(name, variable.datatype, variable.dimensions)
#             dst[name][:] = src[name][:]
#             # copy variable attributes all at once via dictionary
#             dst[name].setncatts(src[name].__dict__)





if __name__ == '__main__':


    lat = np.arange(90.0, -90.0, -1)
    lon = np.arange(180.0, -180.0, -1)

    x, y = np.meshgrid(lon, lat)
    #
    # # 先创建关联的维度，如果不创建关联维度，则创建数据集时，默认赋值range(shape)
    # # writenc('test.nc', 'lat', lat, overwrite=1)
    # # writenc('test.nc', 'lon', lon, overwrite=0)
    #
    lb_writenc('test.nc', 'x', x, dimension=('lat', 'lon'), overwrite=1, dictsdsinfo={'b':1})
    lb_writenc_sdsinfo('test.nc', 'x',  overwrite=0, dictsdsinfo={'bcc':1})
    lb_writenc_fileinfo('test.nc', {'temp':12}, overwrite=0)

    data, info = lb_readnc('test.nc', 'x', dictsdsinfo={})
    print(data, info)
    fileinfo = lb_readnc_fileinfo('test.nc')
    print(fileinfo)
    sdsinfo = lb_readnc_sdsinfo('test.nc', 'x')
    print(sdsinfo)

    # lb_writenc('test.nc', 'y', y, dimension=('lat', 'lon'), overwrite=0,
    #         dictsdsinfo={'b':1})

    # data, fileinfo = lb_readnc('test.nc', 'y', dictfileinfo={})
    # print(data)
    # print(fileinfo)






