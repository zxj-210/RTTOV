#coding:utf-8

import os
import h5py

def lb_readhdf(filename, sdsname, dictsdsinfo = None, dictfileinfo = None) :
    '''
    读取hdf数据集和属性
    :param filename:
    :param sdsname:
    :param dictsdsinfo:
    :param dictfileinfo:
    :return:
    '''
    if not os.path.isfile(filename) :
        print('%s is not exist, will return None' %(filename))
        return None
    else:
        fp = h5py.File(filename, 'r')
        dsetid = fp[sdsname]

        data = dsetid[:]
        if not dictsdsinfo is None :
            for key in dsetid.attrs :
                # print(key)
                dictsdsinfo.update({key : dsetid.attrs[key]})
            fp.close()
            return data, dictsdsinfo

        fp.close()

        return data


def lb_writehdf(filename, sdsname, data, overwrite=True, dictsdsinfo = None, dictfileinfo = None, compression = 9):
    '''
    创建hdf5文件
    :param filename:
    :param sdsname:
    :param data:
    :param overwrite:
    :param dictsdsinfo:
    :param dictfileinfo:
    :param compression:
    :return:
    '''

    try:

        if overwrite :
            fp = h5py.File(filename, 'w')
        else:
            fp = h5py.File(filename, 'r+')
        if not dictfileinfo is None :
            for key in dictfileinfo :
                fp.attrs[key] = dictfileinfo[key]

        dsetid = fp.create_dataset(sdsname, data=data, compression = compression)

        if not dictsdsinfo is None :
            for key in dictsdsinfo :
                dsetid.attrs[key] = dictsdsinfo[key]

        fp.close()
    except BaseException as e :
        print(e)
        return False

    return True


