'''
Created on Nov 20, 2015
updated 11/23/2015

@author: rhealy
'''
import numpy as np

class QuatHeader:

    def __init__(self, **kwargs):
        self.properties = kwargs

    def get_properties(self):
        return self.properties

    def get_property(self, key):
        return self.properties.get(key, None)

def checkData(indata,X,Y,Z,scl,lower,upper):
    badResult = False
    rsltx = (indata[X]*scl)**2
    rslty = (indata[Y]*scl)**2
    rsltz = (indata[Z]*scl)**2
    sumsqrs = np.sqrt(rsltx+rslty+rsltz)
    idx = (sumsqrs > upper) | (sumsqrs < lower)
    if idx.any() :
        badResult = True
    return badResult

def checkQData(indata):
    badResult = False
    rsltx = (indata['ISSQX'])**2
    rslty = (indata['ISSQY'])**2
    rsltz = (indata['ISSQZ'])**2
    rslts = (indata['ISSQS'])**2
    sumsqrs = np.sqrt(rsltx+rslty+rsltz+rslts) - 1.0
    idx = (np.abs(sumsqrs) > 0.01)
    if idx.any() :
        badResult = True
    return badResult

def read_pos_vel_quat(filename):
    import pandas as pd
    import sys
    re=6378137e0    # WGS-84
    lower=re + 2.5e5
    upper=re + 5.e5 # upper and lower bounds for ISS in m wrt WGS-84
    header = {}
    quat = QuatHeader()
    lcnt= 0
    cnt = 0
    for theLine in open(filename,'r') :
        cnt += 1
        #print(lcnt)
        fields = theLine.split(',')
        if len(fields) > 1 and len(fields) < 10 and  'This field is currently' not in fields and 'SecondsSinceEpoch' not in fields:
            header[fields[0]] = ''.join(fields[1]).strip()
            setattr(quat,fields[0],header[fields[0]])
#            print(header[fields[0]])
        if 'SecondsSinceEpoch' in fields:
            lcnt = cnt - 1
    pvq_data = pd.read_csv(filename,skiprows=lcnt,skipinitialspace=True)
    if np.isnan(pvq_data['ISSPOSX']).any() or np.isnan(pvq_data['ISSPOSY']).any() or np.isnan(pvq_data['ISSPOSZ']).any():
        print('Returning; NaN detected in input position')
        print('PLEASE CHECK FILE: {}'.format(filename))
        print('And its source files: \n{}\n{}\n'.format(header['Source CSV file'],header['Subset CSV file']))
        sys.exit()
    if np.isnan(pvq_data['ISSQX']).any() or np.isnan(pvq_data['ISSQY']).any() or np.isnan(pvq_data['ISSQZ']).any():
        print('Returning; NaN detected in input quaternion')
        print('PLEASE CHECK FILE: {}'.format(filename))
        print('And its source files: \n{}\n{}\n'.format(header['Source CSV file'],header['Subset CSV file']))
        sys.exit()
    if checkData(pvq_data,'ISSPOSX','ISSPOSY','ISSPOSZ',.3048,lower,upper):
        print('ERROR: POSITIONS ARE BELOW 200 km above WGS-84 or above 500 km above WGS-84!!!!')
        print('PLEASE CHECK FILE: {}'.format(filename))
        print('And its source files: \n{}\n{}\n'.format(header['Source CSV file'],header['Subset CSV file']))
        sys.exit()
    if checkData(pvq_data,'ISSVELX','ISSVELY','ISSVELZ',.3048,6500.0,9000.0):
        print('ERROR: ISS velocity out of range!')
        print('PLEASE CHECK FILE: {}'.format(filename))
        print('And its source files: \n{}\n{}\n'.format(header['Source CSV file'],header['Subset CSV file']))
        sys.exit()
    if checkQData(pvq_data):
        print('ERROR: Interpolated ISS USGNC quaternions are far from normalized!')
        print('PLEASE CHECK FILE: {}'.format(filename))
        print('And its source files: \n{}\n{}\n'.format(header['Source CSV file'],header['Subset CSV file']))
        sys.exit()
    # Now, the times in the HICO files are in seconds since 2009-Jan01-00:00:00 UTC.
    # N.B. The original files have GPS times in them, and the file interpolate_fields_to_hico_times.pro
    # calculates the GPS second for the above UTC date, and subtracts it form the GPS times in the files.
    # Now, I've verified that the times in the original CSV files are GPS. The YYYYMMDD hhmmss in the files
    # are not UTC, but are instead just GPS time converted. So they differ from UTC by the number of leap
    # seconds there are. So, while my EPOCH date is UTC, NO LEAP SECONDS are included since the epoch. Which means
    # that after 2012/06/30 23:59:59, I have an added second (leap second) that I will need to account for.
    # So, I need to properly count seconds and properly account for leap seconds in order to really get/stay in UTC.
#    print ('Epoch = {} '.format(header['Epoch']))
    if header['Epoch'] != '2009 Jan 01 00:00:00 UTC':
        print('Error: Epoch = {}'.format(header['Epoch']))
        print('In File: {}'.format(filename))
        print('Epoch is not the expected epoch.')
        sys.exit()

    return header,quat,pvq_data

if __name__ == '__main__':

    (quat_info, quat, pvq_data) = read_pos_vel_quat('iss.2013067.0308.063527.L0.12933.20130308205743.hico_pos_vel_quat.csv')
    print('The CSV file is...' + quat_info['Source CSV file'])
    print('Theta = ' + quat_info['Theta (degrees from stowed position)'])
# print('The CSV file is...' + getattr(quat_info,'Source CSV file'))
    print('Distance Unit =' + getattr(quat, 'Distance Unit'))
