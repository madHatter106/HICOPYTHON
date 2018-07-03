#!/usr/bin/env python3
'''
Created on Nov 19, 2015

@author: rhealy
'''
from netCDF4 import Dataset

def proc_hico(ifile,earth_orient_file,leap_sec_file,boresight,nc_grp,calib_grp):
    from hico.hicoLonLatHdr import hicoLonLatHdr
    from hico.astreduc import frac2DegMinSec
    import hico.bore_sight as bore
    import hico.read_pos_vel_quat as rquat
    import numpy as np
    import os,sys
    import struct
    import hico.astreduc as astreduc

    lonHemi = {1:'E',0:'W'}
    latHemi = {1:'N',0:'S'}
    direction = {1:'ascending',0:'descending'}
    byteorder = {'little':0, 'big':1}

    rad2Deg = 180.0/np.pi
    deg2Rad = np.pi/180.0
    theta_stowed=-60.

    try:
        ocvarroot = os.environ['OCVARROOT']
    except KeyError:
        print("OCSSW environement variable OCVARROOT not set.")
        print("Typically it is set to $OCSSW/run/var.")
        sys.exit()


    bs = np.zeros((3))
    bs[0] = boresight[0]*deg2Rad
    bs[1] = boresight[1]*deg2Rad
    bs[2] = boresight[2]*deg2Rad
    r_hico_bs = bore.bore_sight(bs)

    (quat_info, quat, pvq_data) = rquat.read_pos_vel_quat(ifile)
    nsamples = np.int(quat_info['Requested number of pixels'])
#    print('The CSV file is...' + quat_info['Source CSV file'])
#    print('Theta = ' + quat_info['Theta (degrees from stowed position)'])
#    print('The CSV file is...' + getattr(quat_info,'Source CSV file'))
#    print('Distance Unit =' + getattr(quat, 'Distance Unit'))

    filein = "{}/hico/{}".format(ocvarroot,'nutation.dat')
#    print(filein)
    nut_data = astreduc.initReduc(filein)
    astrodate = astreduc.HicoSec2Date(pvq_data['SecondsSinceEpoch'][0])

    yyyy = astrodate.yyyy % 100

    try:
        pmx,pmy,dut1,loda = astreduc.getEDFData(earth_orient_file,yyyy,astrodate.mm,astrodate.dd)
#        print(pmx,pmy,dut1,loda)
    except ValueError:
        print("Earth Data File does not have year={},month={},day={}".format(astrodate.yyyy,astrodate.mm,astrodate.dd))
        sys.exit()

    try:
#        print("Date={}/{}/{} {}:{}:{}".format(astrodate.yyyy,astrodate.mm,astrodate.dd,astrodate.hour,astrodate.min,astrodate.sec))
        lsdat = astreduc.getLSData(leap_sec_file,astrodate.yyyy,astrodate.mm,astrodate.dd)
#        print('lsdat=' + str(lsdat))
    except ValueError:
        print("Leap Sec File does not have year={},month={},day={}".format(astrodate.yyyy,astrodate.mm,astrodate.dd))
        sys.exit()

# Define the rotation from HICO frame of reference to ISS frame of reference
# NOTE: This assumes perfect alignment.

    #r_hico_to_iss = np.zeros((3,3))
    #r_hico_to_iss[0][0] = -1.0
    #r_hico_to_iss[1][1] = -1.0
    #r_hico_to_iss[2][2] =  1.0
    r_hico_to_iss = np.diagflat([-1]*2+[1])
# include bore sight offsets

    r_hico_to_iss = np.dot(r_hico_to_iss,r_hico_bs)
#!!!! BEGIN HICO pointing angles:
# I believe the 1/2 FOV is 3.461515 degrees (from Mike C?)
# A document by Bob Lucke states that the relationship is as follows;
# this was in MJM's HICO email folder in an email from Bill Snyder.
# UNITS ARE DEGREES

    if quat_info['ISS orientation'] == '+XVV':
        imult = 1
    else:
        imult = -1
    theta = np.zeros(nsamples)

    for i in range(0,nsamples):
        frac = imult*(i-255.5)/255.5
        theta[i] = ( -3.487  + 0.035 * frac**2)*frac

    theta = (theta + theta_stowed + float(quat_info['Theta (degrees from stowed position)']))*deg2Rad
# Now theta is in radians and accounts for angle from stowed position
# Now I need to convert to pointing angles in HICO's frame of reference.
# Note that HICO is in the YZ plane of HICO. The stowed angle (+ is right
# hand rule about X, measured from +Z) is in the +Y+Z quadrant, even
# for negative values of  theta...

    svec = np.zeros((3,nsamples))

    svec[0,:] = 0.0
    svec[1,:] = -np.sin(theta[:])
    svec[2,:] =  np.cos(theta[:])

    pm = astreduc.polarm(pmx,pmy)

    #print(pm)

    r_iss = np.array((3),dtype=float)
    q_iss = np.array((4),dtype=float)
    fout = quat_info['Expected Lon/Lat/View angle filename']
    fheader = fout.split('.bil')[0] + '.hdr'
    fbil=open(fout,"wb")
    lines = pvq_data['SecondsSinceEpoch'].size

    #Setup the HICO header file using the hicoLonLatHdr class
    hicoHeader = hicoLonLatHdr(filename=fheader, \
                    samples=nsamples, \
                    lines=lines, \
                    description='HICO geometry file, with calculated positions, and solar & view geometry, on the ground', \
                    sensor_type='HICO-ISS' \
                )
    hicoHeader.set_variable('data type',5)
    hicoHeader.set_variable('interleave','bil')
    hicoHeader.set_variable('wavelength units','{ degrees, degrees, degrees, degrees, degrees, degrees }')
    hicoHeader.set_variable('band names','{ longitude, latitude, view zenith, view azimuth, sol zenith, sol azimuth }')
    hicoHeader.set_variable('byte order',byteorder[sys.byteorder])
    hicoHeader.set_variable('geometry_deltaUT1_s',dut1)
    hicoHeader.set_variable('geometry_deltaT_s',lsdat)
    hicoHeader.set_variable('geometry_xp_rad',pmx)
    hicoHeader.set_variable('geometry_yp_rad',pmy)
    hicoHeader.set_variable('geometry_LOD',loda)
    hicoHeader.set_variable('geometry_bore_sight_params','{ ' + ', '.join( str(v) for v in bs) + ' }')

    for i in range(0,lines):
#    for i in range(0,1):

        r_iss = np.multiply([pvq_data['ISSPOSX'][i],pvq_data['ISSPOSY'][i],pvq_data['ISSPOSZ'][i] ],0.3048e0) # convert to meeters
        q_iss = [ pvq_data['ISSQS'][i],pvq_data['ISSQX'][i],pvq_data['ISSQY'][i],pvq_data['ISSQZ'][i] ]
        rot_body = astreduc.quat_to_rot(q_iss)
        astrodate = astreduc.HicoSec2Date(pvq_data['SecondsSinceEpoch'][i])
        hico_jdate = astreduc.UT2time(astrodate.yyyy,astrodate.mm,astrodate.dd,astrodate.hour,astrodate.min,astrodate.sec,dut1,lsdat)
        #print(hico_jdate.jdut1,hico_jdate.ttdb)
        prec = astreduc.precession(hico_jdate.ttdb)
        DeltaPsi, TrueEps, MeanEps, Omega, nut = astreduc.nutation(hico_jdate.ttdb,nut_data)
        st,stdot,OmegaEarth = astreduc.sidereal(hico_jdate.jdut1,DeltaPsi,TrueEps,Omega,loda,2)

        t_eci_to_ecef=np.dot(pm,np.dot(st,np.dot(nut,prec)))
        # for attitudes from HICO frame to ECEF
        t_hico_to_ecef=np.dot(t_eci_to_ecef,np.dot(rot_body,r_hico_to_iss))
        # The position does not depend on the attitude transformation
        r_ecef=np.dot(t_eci_to_ecef,r_iss)
        # Here, v_ecef are the 512 pointing vectors.
        # v_ecef=matmul(t_hico_to_ecef,v_hico) ! v_hico=[3,NS]
        v_ecef=np.dot(t_hico_to_ecef,svec)
        llh,view_zen,view_az = astreduc.wgs84_intercept(r_ecef,v_ecef)

        sol_zen_az=rad2Deg*astreduc.solar_geometry(astrodate.yyyy,astrodate.mm,astrodate.dd,astrodate.hour,astrodate.min,astrodate.sec,llh[1,:],-llh[0,:])

        lon_out = np.array(llh[0,:])
        lat_out = np.array(llh[1,:])
        viewzen = np.array(view_zen,dtype='f8')
        viewaz  = np.array((360.0 +view_az) % 360.0,dtype='f8')
        solarzen = np.array(sol_zen_az[0,:],dtype='f8')
        solaraz  = np.array(sol_zen_az[1,:],dtype='f8')
        nc_grp['sensor_zenith'][i] = viewzen
        nc_grp['sensor_azimuth'][i] = viewaz
        nc_grp['solar_zenith'][i] = solarzen
        nc_grp['solar_azimuth'][i] = solaraz
        nc_grp['longitudes'][i] = lon_out
        nc_grp['latitudes'][i]  = lat_out
        bilout = np.concatenate((lon_out,lat_out))
        bilout = np.concatenate((bilout,viewzen))
        bilout = np.concatenate((bilout,viewaz))
        bilout = np.concatenate((bilout,solarzen))
        bilout = np.concatenate((bilout,solaraz ))
        myfmt='d'*len(bilout)
        binout = struct.pack(myfmt,*bilout)
        fbil.write(binout)

        # Below are some items for output header, grabbed near center of line
        # Populate the HICO header for output to the header file
        if i == lines/2-1:
            hicoHeader.set_variable('date', '{ ' + str(astrodate.yyyy)+','+str(astrodate.mm)+','+str(astrodate.dd) + ' }')
            hicoHeader.set_variable('time', '{ ' + str(astrodate.hour)+','+str(astrodate.min)+','+str(astrodate.sec) + ' }')
            hicoHeader.set_variable('image_center_long','{ ' + ', '.join( str(v) for v in frac2DegMinSec(np.abs(lon_out[254]),0.,0.)) + ' }')
            hicoHeader.set_variable('image_center_lat', '{ ' + ', '.join( str(v) for v in frac2DegMinSec(np.abs(lat_out[254]),0.,0.)) + ' }')
            hicoHeader.set_variable('image_center_long_hem', lonHemi[(lon_out[254]>=0)])
            hicoHeader.set_variable('image_center_lat_hem', latHemi[(lat_out[254]>=0)])
            hicoHeader.set_variable('image_center_zenith_ang','{ ' + ', '.join( str(v) for v in frac2DegMinSec((viewzen[254]+viewzen[255])/2,0.,0.)) + ' }')
            hicoHeader.set_variable('image_center_azimuth_ang','{ ' + ', '.join( str(v) for v in frac2DegMinSec((viewaz[254]+viewaz[255])/2,0.,0.)) + ' }')
            hicoHeader.set_variable('geometry_ISS_Z_direction',direction[(pvq_data['ISSVELZ'][i]>=0)]) #! a bit crude, based on ECI velocity
            r_iss_new = np.tile(r_iss[np.newaxis],1).T
            llh = astreduc.ecef2latlon(r_iss_new)
            hicoHeader.set_variable('sensor_altitude',llh[2,0]/1000.)

    fbil.close()

    # Create the history for output into the header file using the information from the
    # pos_vel_quat input csv file and then write out the header file.

    hicoHeader.createHistory(quat_info)
    hicoHeader.writeHicoENVIHeaderFile()

    calib_grp.hico_orientation_from_quaternion = quat_info['ISS orientation']

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description='''\
      Generate a HICO geometry file.

            ''', add_help=True)
    parser.add_argument('-ifile', nargs=1, type=str, help=' iss*pos_vel_quat.csv input file (Must be a CSV file) ')
    parser.add_argument('-lfile', nargs=1, type=str, help=' leapsec file ')
    parser.add_argument('-efile', nargs=1, type=str, help=' earth-orient file ')
    parser.add_argument('-boresight', nargs=3, type=float, default=([0, 0, 0]), help=('Process Bore Sight Parameters'))

    args = parser.parse_args('-ifile iss.2013067.0308.063527.L0.12933.20130308205743.hico_pos_vel_quat.csv -efile /home/rhealy/ocsswn/run/var/hico/finals.data -lfile /home/rhealy/ocsswn/run/var/modis/leapsec.dat -boresight -0.9957 0.0268 -0.0128'.split())

    ifile = ''.join(args.ifile)
    earth_orient_file = ''.join(args.efile)
    leap_sec_file     = ''.join(args.lfile)
    boresight        = args.boresight
    rootgrp = Dataset('test_hico.nc','w')
    rootgrp.createDimension('scan_lines', 2000)
    rootgrp.createDimension('samples',  512)
    #
    # L2gen is expecting these variables under group "navigation" for l2gen
    #
    navigation = rootgrp.createGroup("navigation")

    senz = navigation.createVariable('sensor_zenith','f4',('scan_lines','samples',))
    solz = navigation.createVariable('solar_zenith','f4',('scan_lines','samples',))
    sena = navigation.createVariable('sensor_azimuth','f4',('scan_lines','samples',))
    sola = navigation.createVariable('solar_azimuth','f4',('scan_lines','samples',))
    lon  = navigation.createVariable('longitudes','f4',('scan_lines','samples',))
    lat  = navigation.createVariable('latitudes','f4',('scan_lines','samples',))
    senz.units = 'degrees'
    senz.valid_min = -180
    senz.valid_max =  180
    senz.long_name = 'sensor zenith'
    solz.units = 'degrees'
    solz.valid_min = -180
    solz.long_name = 'solar zenith'
    solz.valid_max =  180
    sena.units = 'degrees'
    sena.valid_min = -180
    sena.valid_max =  180
    sena.long_name = 'sensor azimuth'
    sola.units = 'degrees'
    sola.valid_min = -180
    sola.valid_max =  180
    sola.long_name = 'solar azimuth'
    lon.units = 'degrees'
    lon.valid_min = -180
    lon.valid_max =  180
    lon.long_name = 'longitude'
    lat.units = 'degrees'
    lat.valid_min = -180
    lat.valid_max =  180
    lat.long_name = 'latitude'
    # orientation needs to be set in group = "metadata/HICO/Calibration for l2gen"
    metadata        = rootgrp.createGroup("metadata")
    hico            = metadata.createGroup("HICO")
    calibration     = hico.createGroup("Calibration")
    proc_hico(ifile,earth_orient_file,leap_sec_file,boresight,navigation,calibration)

    rootgrp.close()
