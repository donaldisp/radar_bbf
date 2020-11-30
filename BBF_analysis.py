# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:06:52 2020

@author: Donaldi Permana
"""

import wradlib as wrl
import matplotlib.pyplot as pl
import matplotlib as mpl
import warnings
warnings.filterwarnings('ignore')
# try:
#     get_ipython().magic("matplotlib inline")
# except:
#     pl.ion()
import numpy as np
from osgeo import osr
import datetime as dt
# import gc
import os

def process_radar_data_gematronik(filename):
    print ('Processing Gematronik Radar Data :'+filename)
    
    rmax = 0
    degtometer = 111229.
    #filename = 'I:\Radar_Integrasi\BKS20170919210000dBZ.vol'
    raw = wrl.io.rainbow.read_rainbow(filename)
    # newdir = filename+'.dir'
    # try:
    #     os.mkdir(newdir, 755 );
    # except:
    #     print (newdir+' already exists')
    # this is the radar position tuple (longitude, latitude, altitude)
    try:
        llon = float(raw['volume']['sensorinfo']['lon'])
        llat = float(raw['volume']['sensorinfo']['lat'])
        lalt = float(raw['volume']['sensorinfo']['alt'])
        lbeamwidth = float(raw['volume']['sensorinfo']['beamwidth'])
        lidname = raw['volume']['sensorinfo']['@id']
        lsensortype = raw['volume']['sensorinfo']['@type']
        lsensorname = raw['volume']['sensorinfo']['@name']
    except: #aceh
        llon = float(raw['volume']['radarinfo']['@lon'])
        llat = float(raw['volume']['radarinfo']['@lat'])
        lalt = float(raw['volume']['radarinfo']['@alt'])
        lbeamwidth = float(raw['volume']['radarinfo']['beamwidth'])
        lidname = raw['volume']['radarinfo']['@id']
        lsensortype = ''
        lsensorname = raw['volume']['radarinfo']['name']
        
    sitecoords = (llon, llat,lalt)
    beamwidth = lbeamwidth
    sitename = lsensorname + ' ' + lidname
    
    # dtime = raw['volume']['@datetime']
    # strdt = dt.datetime.strptime(dtime.replace('T',''), "%Y-%m-%d%H:%M:%S")

    # containers to hold Cartesian bin coordinates and data
    xyz, data = np.array([]).reshape((-1, 3)), np.array([])
    
    nelevangle = np.size(raw['volume']['scan']['slice'])
    elevation = np.array([])
    nrays = np.array([])
    nbins = np.array([])
    range_res = np.array([])
    
    # iterate over 14 elevation angles
    for i in range(nelevangle):
    #    try:    
            try:
                angle = round(float(raw['volume']['scan']['slice'][i]['posangle']),2)
                if angle in elevation:
                    continue
                else:                        
                    elevation = np.append(elevation,angle) # in degree
            except:                    
                continue
            
            print ('Reading and plotting vertical scanning angle '+str(angle)+' deg ...')
                            
            # get azimuthal data
            try: #aceh
                azi = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo']['data']
                azidepth = float(raw['volume']['scan']['slice'][i]
                                 ['slicedata']['rayinfo']['@depth'])
                azirange = float(raw['volume']['scan']['slice'][i]
                                 ['slicedata']['rayinfo']['@rays'])
                
            except:
                azi0 = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][0]['data']
                azi1 = raw['volume']['scan']['slice'][i]['slicedata']['rayinfo'][1]['data']
                azi = (azi0/2) + (azi1/2)
                del azi0, azi1
                
                azidepth = float(raw['volume']['scan']['slice'][i]
                                 ['slicedata']['rayinfo'][0]['@depth'])
                azirange = float(raw['volume']['scan']['slice'][i]
                                 ['slicedata']['rayinfo'][0]['@rays'])
            try:
                azires = float(raw['volume']['scan']['slice'][i]['anglestep'])
            except:
                azires = float(raw['volume']['scan']['slice'][0]['anglestep'])
         
            azi = (azi * azirange / 2**azidepth) * azires
            
            nrays = np.append(nrays,len(azi))
            
            flag=0
            if np.size(azi) >= 999:
                flag=2
                azi = azi/3
                for ii in range(int(np.floor(np.size(azi)/3))):
                    azi[ii] = azi[3*ii]+azi[3*ii+1]+azi[3*ii+2]
                azi = azi[range(int(np.floor(np.size(azi)/3)))]
            elif np.size(azi) >= 500:
                flag=1
                azi = azi/2
                for ii in range(int(np.floor(np.size(azi)/2))):
                    azi[ii] = azi[2*ii]+azi[2*ii+1]
                azi = azi[range(int(np.floor(np.size(azi)/2)))]                
            
            # create range array
            try:
                stoprange = float(raw['volume']['scan']['slice'][i]['stoprange'])
                rangestep = float(raw['volume']['scan']['slice'][i]['rangestep'])
            except:
                stoprange = float(raw['volume']['scan']['slice'][0]['stoprange'])
                rangestep = float(raw['volume']['scan']['slice'][0]['rangestep'])
                            
#                if stoprange_max < stoprange:
#                stoprange_max = stoprange
            
            r = np.arange(0, stoprange, rangestep)*1000

            nbins = np.append(nbins,len(r))
            
            range_res = np.append(range_res,rangestep*1000)
            
#             if i==0:
#                 rr=r
                
#             if rmax < r.max():
#                 rmax = r.max()

#             # get reflectivity data
#             data_ = raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['data']
#             datadepth = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@depth'])
#             datamin = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@min'])
#             datamax = float(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@max'])
#             data_ = datamin + data_ * (datamax - datamin) / 2 ** datadepth
            
#             if flag==2:
#                 data_ = data_/3
#                 for jj in range(int(np.floor(np.size(data_[:,1])/3))):
#                     data_[jj,:] = data_[3*jj,:] + data_[3*jj+1,:] + data_[3*jj+2,:]
#                 data_ = data_[range(int(np.floor(np.size(data_[:,1])/3))),:]                
#             elif flag==1:
#                 data_ = data_/2
#                 for jj in range(int(np.floor(np.size(data_[:,1])/2))):
#                     data_[jj,:] = data_[2*jj,:] + data_[2*jj+1,:]
#                 data_ = data_[range(int(np.floor(np.size(data_[:,1])/2))),:]
            
#             # get annotation data
#             unit = raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@type']
#             time = raw['volume']['scan']['slice'][i]['slicedata']['@time']
#             date = raw['volume']['scan']['slice'][i]['slicedata']['@date']
                        
# #            nrays = int(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@rays']) # rays
# #            nbins = int(raw['volume']['scan']['slice'][i]['slicedata']['rawdata']['@bins'])    
            
            # If len(azi) == 447 will generate error in wrl.ipol.interpolate_polar
            # "ValueError: operands could not be broadcast together with shapes"
            if len(azi) == 175:
                azi = azi[:-1]
#                 data_ = data_[:-1,:]
                
#             delta = len(r) - len(np.transpose(data_))
#             if delta > 0:
#                 r = r[:-delta]
                
#             clutter = wrl.clutter.filter_gabella(data_, wsize=5, thrsnorain=0.,tr1=6., n_p=8, tr2=1.3)
#             data_no_clutter = wrl.ipol.interpolate_polar(data_, clutter)
#             del data_, clutter
            
# #                pia = wrl.atten.correctAttenuationHJ(data_no_clutter, a_max = 4.565e-5, b = 0.73125, n=1, mode = 'cap', thrs_dBZ = 100.0, max_PIA = 4.82)
                      
# # =============================================================================
# #                 pia = wrl.atten.correct_attenuation_constrained(data_no_clutter, a_max=1.67e-4,
# #                                             a_min=2.33e-5, n_a=100,
# #                                             b_max=0.7, b_min=0.65,
# #                                             n_b=6, gate_length=1.,
# #                                             constraints=
# #                                             [wrl.atten.constraint_dbz,wrl.atten.constraint_pia],
# #                                             constraint_args=
# #                                             [[59.0],[20.0]])
# # 
# #                 data_no_clutter1 = data_no_clutter + pia
# # =============================================================================
            
#             pia = wrl.atten.correct_attenuation_hb(data_no_clutter,coefficients = dict(a=4.57e-5, b=0.731, gate_length=1.0),mode="warn",thrs=59.)
#             pia[pia > 4.8] = 4.8
                            
#             data_attcorr = data_no_clutter + pia
#             #data_attcorr = data_no_clutter
#             del data_no_clutter, pia
            
#             #azimuths = np.arange(0, 360., 360./nrays) # in degrees
#             #ranges = np.arange(0, nbins*rscale, rscale) # in meters
          
#             polargrid = np.meshgrid(r, azi)
# #                lon, lat, alt = wrl.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], angle, sitecoords)
            
#             coords = wrl.georef.polar.spherical_to_proj(polargrid[0], polargrid[1], elevation[i], sitecoords)
#             lon = coords[:,:,0]
#             lat = coords[:,:,1]
#             alt = coords[:,:,2]
#             gk3 = wrl.georef.epsg_to_osr(4326)
#             ae = wrl.georef.create_osr("aeqd", lon_0=llon, lat_0=llat)
#             x_, y_ = wrl.georef.reproject(lon, lat, projection_target=gk3)
#             x, y = wrl.georef.reproject(lon, lat, projection_target=ae)
#             del polargrid, lon, lat, alt
                           
            
#             if i==0:
#                 # define horizontal resolution in this part
#                 res = 500. # meters
#                 ngrid = int(np.floor(((x_.max()-x_.min())/(res/degtometer))+1))
#                 xgrid_ = np.linspace(x_.min(), x_.max(), ngrid)
#                 ygrid_ = np.linspace(y_.min(), y_.max(), ngrid)                
#                 xgrid = np.linspace(x.min(), x.max(), ngrid)
#                 ygrid = np.linspace(y.min(), y.max(), ngrid)
                
#                 xgridmin = xgrid_.min()
#                 xgridmax = xgrid_.max()
#                 ygridmin = ygrid_.min()
#                 ygridmax = ygrid_.max() 
#                 scan_data = np.array([]).reshape((-1, len(ygrid)))
            
#             grid_xy_ = np.meshgrid(xgrid_, ygrid_)                   
#             grid_xy = np.meshgrid(xgrid, ygrid)
#             xx = grid_xy_[0]
#             yy = grid_xy_[1]
#             grid_xy = np.vstack((grid_xy[0].ravel(), grid_xy[1].ravel())).transpose()
            
# #                print 'lon '+str(xgridmin)+ ' to '+str(xgridmax)
# #                print 'lat '+str(ygridmin)+ ' to '+str(ygridmax)
            
#             xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
#             gridded = wrl.comp.togrid(xy, grid_xy, r[np.size(r)-1], np.array([x.mean(), y.mean()]), data_attcorr.ravel(), wrl.ipol.Nearest)
#             gridded = np.ma.masked_invalid(gridded).reshape((len(xgrid), len(ygrid)))
            
#             scan_data = np.vstack((scan_data,gridded))

# #                print 'len '+str(len(gridded))

#             del x, y, grid_xy, xy
#             # xyz_ = wrl.vpr.volcoords_from_polar(sitecoords, angle, azi, r, proj)
#             # get the scan data for this elevation
#             #   here, you can do all the processing on the 2-D polar level
#             #   e.g. clutter elimination, attenuation correction, ...
#         #    data_ = what["offset"] + what["gain"] * raw[
#         #        "dataset%d/data1/data" % (i + 1)]
#             # transfer to containers
#         #    xyz, data = np.vstack((xyz, xyz_)), np.append(data, data_.ravel())
        
#             # xyz, data = np.vstack((xyz, xyz_)), np.append(data, data_attcorr.ravel())
#             # del xyz_, data_attcorr
#     #    except:
#     #        print("\nError data in vertical elevation scanning - %s deg" % elevation)
#     #        print('\nError data in vertical elevation scanning - ' + str(elevation) + ' deg')
#             gc.collect();
#             del gc.garbage[:]
            
    return sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename 

def process_radar_data_eec(filename):
    print ('Processing EEC Radar Data :'+filename)
    
    rmax = 0
    degtometer = 111229.
    #filename = 'I:\Radar_Integrasi\\2020\EEC\JAK-20200501-140214-PPIVol.nc'
    raw = wrl.io.read_generic_netcdf(filename)
#    newdir = filename+'.dir'
#    try:
#        os.mkdir(newdir, 0755 );
#    except:
#        print newdir+' already exists'
    # this is the radar position tuple (longitude, latitude, altitude)
    
    llon = float(raw['variables']['longitude']['data'])
    llat = float(raw['variables']['latitude']['data'])
    lalt = float(raw['variables']['altitude']['data'])
    
    # radar center coordinate correction
    if 'JAK' in filename:
       llat = -6.1715
       llon = 106.6466

    lbeamwidth = float(raw['variables']['radar_beam_width_h']['data'])
    lidname = raw['site_name']
    lsensorname = raw['institution']
        
    sitecoords = (llon, llat,lalt)
    beamwidth = lbeamwidth
    sitename = lsensorname + ' ' + lidname

   # containers to hold Cartesian bin coordinates and data
    xyz, data = np.array([]).reshape((-1, 3)), np.array([])
    
    dtime = str(raw['variables']['time_coverage_start']['data'])
    dtime = dtime.replace('Z','')
    # strdt = dt.datetime.strptime(dtime.replace('T',''), "%Y-%m-%d%H:%M:%S")
    
    nelevangle = np.size(raw['variables']['fixed_angle']['data'])
    # datetime = str(raw['variables']['time_coverage_start']['data'])
    
    sweep_start_idx = raw['variables']['sweep_start_ray_index']['data']
    sweep_end_idx = raw['variables']['sweep_end_ray_index']['data']

    # create range array
    r_all = raw['variables']['range']['data']
    rr = r_all
    
    try:        
        if raw['gates_vary'] == 'true':
            ray_n_gates = raw['variables']['ray_n_gates']['data']
            ray_start_index = raw['variables']['ray_start_index']['data']
            flag = 'true'
        else:
            flag = 'false'
    except:
        if raw['n_gates_vary'] == 'true':
            ray_n_gates = raw['variables']['ray_n_gates']['data']
            ray_start_index = raw['variables']['ray_start_index']['data']
            flag = 'true'
        else:
            flag = 'false'
    
    elevation = np.array([])    
    nrays = np.array([])
    nbins = np.array([])
    range_res = rr
    
    # iterate over 9 elevation angles
    for i in range(nelevangle):

            elevation = np.append(elevation,float('{0:.1f}'.format(raw['variables']['fixed_angle']['data'][i]))) # in degree
    
            print ('Reading vertical scanning angle '+str(elevation[i])+' deg ...')
            
            # get azimuthal data        
            azi = raw['variables']['azimuth']['data'][sweep_start_idx[i]:sweep_end_idx[i]]
           
            # # get reflectivity data
            # if flag == 'false':
            #     try:                    
            #         data_ = raw['variables']['DBZH']['data'][sweep_start_idx[i]:sweep_end_idx[i], :]
            #     except:                    
            #         data_ = raw['variables']['UH']['data'][sweep_start_idx[i]:sweep_end_idx[i], :]
            #     try:
            #         datav_ = raw['variables']['VELH']['data'][sweep_start_idx[i]:sweep_end_idx[i], :]
            #     except:
            #         pass
                
            #     # create range array
            #     r = r_all
            # else: #flag = true
            #     data_ = np.array([])
            #     datav_ = np.array([])
            #     n_azi = sweep_end_idx[i]-sweep_start_idx[i]
            #     try:
            #         for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
            #             data_ = np.append(data_,raw['variables']['DBZH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
            #         data_ = data_.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            #     except:
            #         for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
            #             data_ = np.append(data_,raw['variables']['UH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
            #         data_ = data_.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            #     try:
            #         for ll in range(sweep_start_idx[i],sweep_end_idx[i]):
            #             datav_ = np.append(datav_,raw['variables']['VELH']['data'][ray_start_index[ll]:ray_start_index[ll+1]])
            #         datav_ = datav_.reshape((n_azi,ray_n_gates[sweep_start_idx[i]]))
            #     except:
            #         pass

            #     # create range array
            #     r = r_all[0:ray_n_gates[sweep_start_idx[i]]]
            
            rmax = r_all.max()            
            
            # # get annotation data
            # try:
            #     unit = raw['variables']['DBZH']['units']
            # except:
            #     unit = raw['variables']['UH']['units']
            
            # If len(azi) == 447 will generate error in wrl.ipol.interpolate_polar
            # "ValueError: operands could not be broadcast together with shapes"
            if len(azi) == 447:
                azi = azi[:-1]
            #     data_ = data_[:-1,:]                

            # # Correct the reflectivity data with radial velocity data (V)
            # if ('JAK' in filename or 'MED' in filename or 'PEK' in filename or 'TAR' in filename):
            #     if (elevation[i] < 1.5): #only scan angle less than 2 degree are corrected
            #         try:
            #             data_[(datav_==-32768.)]=-32768. 
            #         except:
            #             pass

            nrays = np.append(nrays,len(azi))
            nbins = np.append(nbins,len(rr))
            
            # # range_res = np.append(range_res,rangestep*1000)

            # #clutter = wrl.clutter.filter_gabella(data_, tr1=6, n_p=6, tr2=1.3, rm_nans=False)
            # clutter = wrl.clutter.filter_gabella(data_, tr1=6., n_p=8, tr2=1.3, wsize=5, thrsnorain=0.)
            # #clutter = wrl.clutter.filter_gabella_a(data_, tr1=10., wsize=9)
            # data_no_clutter = wrl.ipol.interpolate_polar(data_, clutter)
            
            # #pia = wrl.atten.correctAttenuationKraemer(data_no_clutter)
            # # 0.9.0  
            # #pia = wrl.atten.correctAttenuationHJ(data_no_clutter, a_max = 4.565e-5, b = 0.73125, n=1, mode = 'cap', thrs_dBZ = 100.0, max_PIA = 4.82)

            # pia = wrl.atten.correct_attenuation_constrained(data_no_clutter, a_max=1.67e-4,
            #                             a_min=2.33e-5, n_a=100,
            #                             b_max=0.7, b_min=0.65,
            #                             n_b=6, gate_length=1.,
            #                             constraints=
            #                             [wrl.atten.constraint_dbz,
            #                             wrl.atten.constraint_pia],
            #                             constraint_args=
            #                             [[59.0],[20.0]])
            # data_attcorr = data_no_clutter + pia
            # data_attcorr[data_attcorr<0] = np.nan
            # #data_attcorr = data_no_clutter
            # del data_no_clutter, pia
            
            # # # derive 3-D Cartesian coordinate tuples
            # # xyz_ = wrl.vpr.volcoords_from_polar(sitecoords, elevation[i], azi, r, proj)
            
            # # # transfer to containers
            # # xyz, data = np.vstack((xyz, xyz_)), np.append(data, data_attcorr.ravel())
           
            # polargrid = np.meshgrid(r, azi)
            # # 0.9.0      
            # #lon, lat, alt = wrl.georef. polar2lonlatalt_n(polargrid[0], polargrid[1], elevation[i], sitecoords)

            # #1.0.0          
            # #lon, lat, alt = wrl.georef.polar.spherical_to_xyz(r, azi, elevation[i], sitecoords)
            # coords = wrl.georef.polar.spherical_to_proj(polargrid[0], polargrid[1], elevation[i], sitecoords)
            # lon = coords[:,:,0]
            # lat = coords[:,:,1]
            # alt = coords[:,:,2]
            # gk3 = wrl.georef.epsg_to_osr(4326)
            # ae = wrl.georef.create_osr("aeqd", lon_0=llon, lat_0=llat)
            # x_, y_ = wrl.georef.reproject(lon, lat, projection_target=gk3)
            # x, y = wrl.georef.reproject(lon, lat, projection_target=ae)
            # del polargrid, lon, lat, alt
          
            # if i==0:
            #      # define horizontal resolution in this part
            #     res = 500. # meters
            #     ngrid = int(np.floor(((x_.max()-x_.min())/(res/degtometer))+1))
            #     xgrid_ = np.linspace(x_.min(), x_.max(), ngrid)
            #     ygrid_ = np.linspace(y_.min(), y_.max(), ngrid)
            #     xgrid = np.linspace(x.min(), x.max(), ngrid)
            #     ygrid = np.linspace(y.min(), y.max(), ngrid)   
                
            #     xgridmin = xgrid_.min()
            #     xgridmax = xgrid_.max()
            #     ygridmin = ygrid_.min()
            #     ygridmax = ygrid_.max()
            #     scan_data = np.array([]).reshape((-1, len(ygrid)))                
            
            # grid_xy_ = np.meshgrid(xgrid_, ygrid_)                   
            # grid_xy = np.meshgrid(xgrid, ygrid)
            # xx = grid_xy_[0]
            # yy = grid_xy_[1]
            # grid_xy = np.vstack((grid_xy[0].ravel(), grid_xy[1].ravel())).transpose()
            
            # xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
            # gridded = wrl.comp.togrid(xy, grid_xy, r[np.size(r)-1], np.array([x.mean(), y.mean()]), data_attcorr.ravel(), wrl.ipol.Nearest)
            # gridded = np.ma.masked_invalid(gridded).reshape((len(xgrid), len(ygrid)))
                 
            # scan_data = np.vstack((scan_data,gridded))
            # #gridded[gridded==0] = np.nan
            
            # del x, y, grid_xy, xy
            # gc.collect();
            # del gc.garbage[:]
            
    return sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename 

def process_radar_data_baron(filename):
    print ('Processing Baron Radar Data :'+filename)

    rmax = 0
    degtometer = 111229.
    #filename = 'I:\Radar_Integrasi\\2020\BARON\KPG20181211041000.vol'
    raw, metaraw = wrl.io.read_gamic_hdf5(filename)
    # global newdir, dtime, strdt
#    newdir = filename+'.dir'
#    try:
#        os.mkdir(newdir, 0755 );
#    except:
#        print newdir+' already exists'
    # this is the radar position tuple (longitude, latitude, altitude)
    llon = float(metaraw['VOL']['Longitude'])
    llat = float(metaraw['VOL']['Latitude'])
    lalt = float(metaraw['VOL']['Height'])
    
    try:
        lbeamwidth = float(metaraw['SCAN0']['pulse_width_us'])
    except:
        lbeamwidth = float(metaraw['SCAN0']['pulse_width'])
    # lidname = raw['site_name']
    # lsensorname = raw['institution']
        
    sitecoords = (llon, llat,lalt)
    beamwidth = lbeamwidth
    # sitename = lsensorname + ' ' + lidname

    sitename = ''
    sitecoords = (llon, llat,lalt)
    
    # define your cartesian reference system
    #proj = wrl.georef.create_osr(32632)
    # proj = osr.SpatialReference()
    
    # # select spatial reference for Indonesia from EPSG (http://www.spatialreference.org/ref/?search=indonesia)
    # if (llat >=0): # north hemisphere
    #     if (llon < 96):
    #         epsg_proj = 23846
    #     elif (llon >= 96 and llon < 102):
    #         epsg_proj = 23847
    #     elif (llon >= 102 and llon < 108):
    #         epsg_proj = 23848
    #     elif (llon >= 108 and llon < 114):
    #         epsg_proj = 23849
    #     elif (llon >= 114 and llon < 120):
    #         epsg_proj = 23850
    #     elif (llon >= 120 and llon < 126):
    #         epsg_proj = 23851
    #     elif (llon >= 126 and llon < 132):
    #         epsg_proj = 23852
    #     else:
    #         epsg_proj = 4238
    # else: # south hemisphere
    #     if (llon >= 96 and llon < 102):
    #         epsg_proj = 23887
    #     elif (llon >= 102 and llon < 108):
    #         epsg_proj = 23888
    #     elif (llon >= 108 and llon < 114):
    #         epsg_proj = 23889
    #     elif (llon >= 114 and llon < 120):
    #         epsg_proj = 23890
    #     elif (llon >= 120 and llon < 126):
    #         epsg_proj = 23891
    #     elif (llon >= 126 and llon < 132):
    #         epsg_proj = 23892
    #     elif (llon >= 132 and llon < 138):
    #         epsg_proj = 23893
    #     elif (llon >= 138):
    #         epsg_proj = 23894 
    #     else:
    #         epsg_proj = 4238
            
    # proj.ImportFromEPSG(epsg_proj) 

#    dtime = metaraw['SCAN0']['Time']
#    try:
#        strdt = dt.datetime.strptime(dtime[:-1].replace('T',''), "%Y-%m-%d%H:%M:%S")
#        dtime = dtime[:-1]
#    except:
#        strdt = dt.datetime.strptime(dtime[:-5].replace('T',''), "%Y-%m-%d%H:%M:%S")
#        dtime = dtime[:-5]

    # readtime from filename - there is a problem with datetime in YOG radar
    dtime = filename[-18:-4]
    strdt = dt.datetime.strptime(dtime, "%Y%m%d%H%M%S")

    # containers to hold Cartesian bin coordinates and data
    xyz, data = np.array([]).reshape((-1, 3)), np.array([])
    
    nelevangle = len(raw)
    elevation = np.array([])
    nrays = np.array([])
    nrays = np.append(nrays,float(metaraw['SCAN0']['ray_count']))
    nbins = np.array([])
    nbins = np.append(nbins,float(metaraw['SCAN0']['bin_count']))
    range_res = metaraw['SCAN0']['r']
    
    rmax = 0
    
    # iterate over 14 elevation angles
    # for i in range(nelevangle):
    #         txtscan = 'SCAN'+str(i)
    
    #         elevation = np.append(elevation,float('{0:.1f}'.format(metaraw[txtscan]['elevation']))) # in degree
            
    #         print 'Reading vertical scanning angle '+str(elevation[i])+' deg ...'
            
    #         # get azimuthal data
    #         azi = metaraw[txtscan]['az']
         
    #         # create range array
    #         r = metaraw[txtscan]['r']
    #         r1 = r
            
    #         if len(r1) > 2000:                
    #             r = r/2
    #             for ii in range(np.size(r)/2):
    #                 r[ii] = r[2*ii]+r[2*ii+1]
    #             r = r[range(np.size(r)/2)]
            
    #         rangestep = metaraw[txtscan]['range_step']
    #         rangestep = r[1]-r[0]
            
    #         if i==0:
    #             rr=r

    #         # get reflectivity data
    #         data_ = raw[txtscan]['Z']['data']
            
    #         if len(r1) > 2000:
    #             data_ = data_/2
    #             for jj in range(int(np.floor(np.size(data_[1,:])/2))):
    #                 data_[:,jj] = data_[:,2*jj] + data_[:,2*jj+1]
    #             data_ = data_[:,range(int(np.floor(np.size(data_[1,:])/2)))]
            
    #         if rmax < r.max():
    #             rmax = r.max()
            
    #         # get annotation data
    #         unit = 'dBZ'                
    #         datetime1 = metaraw[txtscan]['Time']
    #         try:
    #             strdt1 = dt.datetime.strptime(datetime1[:-1].replace('T',''), "%Y-%m-%d%H:%M:%S")
    #             datetime = datetime1[:-1]
    #         except:
    #             strdt1 = dt.datetime.strptime(datetime1[:-5].replace('T',''), "%Y-%m-%d%H:%M:%S")
    #             datetime = datetime1[:-5]
                
    #         # derive 3-D Cartesian coordinate tuples
    #         xyz_ = wrl.vpr.volcoords_from_polar(sitecoords, elevation[i], azi, r, proj)    

    #         clutter = wrl.clutter.filter_gabella(data_, tr1=6, n_p=6, tr2=1.3, rm_nans=False)
    #         data_no_clutter = wrl.ipol.interpolate_polar(data_, clutter)
    #         del data_, clutter
    #         pia = wrl.atten.correctAttenuationHJ(data_no_clutter, a_max = 4.565e-5, b = 0.73125, n=1, mode = 'cap', thrs_dBZ = 100.0, max_PIA = 4.82)
    #         data_attcorr = data_no_clutter + pia
    #         #data_attcorr = data_no_clutter
            
    #         # transfer to containers
    #         xyz, data = np.vstack((xyz, xyz_)), np.append(data, data_attcorr.ravel())            
            
    #         del data_no_clutter, pia
            
    #         polargrid = np.meshgrid(r, azi)
    #         lon, lat, alt = wrl.georef.polar2lonlatalt_n(polargrid[0], polargrid[1], elevation[i], sitecoords)
    #         gk3 = wrl.georef.epsg_to_osr(4326)
    #         ae = wrl.georef.create_osr("aeqd", lon_0=llon, lat_0=llat)
    #         x_, y_ = wrl.georef.reproject(lon, lat, projection_target=gk3)
    #         x, y = wrl.georef.reproject(lon, lat, projection_target=ae)
    #         del polargrid, lon, lat, alt
            
    #         if i==0:
    #              # define horizontal resolution in this part
    #             res = 500. # meters
    #             ngrid = np.floor(((x_.max()-x_.min())/(res/degtometer))+1)
    #             xgrid_ = np.linspace(x_.min(), x_.max(), ngrid)
    #             ygrid_ = np.linspace(y_.min(), y_.max(), ngrid)
    #             xgrid = np.linspace(x.min(), x.max(), ngrid)
    #             ygrid = np.linspace(y.min(), y.max(), ngrid)
                
    #             xgridmin = xgrid_.min()
    #             xgridmax = xgrid_.max()
    #             ygridmin = ygrid_.min()
    #             ygridmax = ygrid_.max() 
    #             scan_data = np.array([]).reshape((-1, len(ygrid)))
                
    #         grid_xy_ = np.meshgrid(xgrid_, ygrid_)
    #         grid_xy = np.meshgrid(xgrid, ygrid)
    #         xx = grid_xy_[0]
    #         yy = grid_xy_[1]
    #         grid_xy = np.vstack((grid_xy[0].ravel(), grid_xy[1].ravel())).transpose()
            
    #         xy=np.concatenate([x.ravel()[:,None],y.ravel()[:,None]], axis=1)
    #         gridded = wrl.comp.togrid(xy, grid_xy, r[np.size(r)-1], np.array([x.mean(), y.mean()]), data_attcorr.ravel(), wrl.ipol.Nearest)
    #         gridded = np.ma.masked_invalid(gridded).reshape((len(xgrid), len(ygrid)))

    #         if i==0:
    #             scan_data_null = np.ma.masked_invalid(gridded).reshape((len(xgrid), len(ygrid)))
    #             scan_data_null[:] = np.nan
            
    #         try:                
    #             scan_data = np.vstack((scan_data,gridded))
    #         except:                    
    #             scan_data = np.vstack((scan_data,scan_data_null))

    #         del x, y, grid_xy, xy
            
    #         gc.collect();
    #         del gc.garbage[:]
            
    # return sitecoords, elevation, xgrid_, ygrid_, scan_data, dtime, rr, rmax, proj, xyz, data

    return sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename

# just a little helper function to style x and y axes of our maps
def annotate_map(ax, cm=None, title=""):
    ticks = (ax.get_xticks()/1000).astype(np.int)
    # ax.set_xticklabels(ticks)
    ax.set_xticklabels('')
    ticks = (ax.get_yticks()/1000).astype(np.int)
    ax.set_yticklabels(ticks)
    # ax.set_xlabel("Kilometers")
    ax.set_ylabel("Kilometers")
    if not cm is None:
        pl.colorbar(cm, ax=ax)
    if not title=="":
        ax.set_title(title, fontdict={'fontsize': 10, 'fontweight': 'medium'})
    ax.grid()

def plotBBFanalysis(sitecoords1, nrays1, nbins1, beamwidth1, range_res1, sitename1, ds, elev1, elev2, elev3):    
        
    # sitecoords = (7.071663, 50.73052, 99.5)
    #sitecoords = (102.341280, -3.858720, 45.0) # Bengkulu
    sitecoords = sitecoords1
    nrays = int(nrays1[0]) # number of rays
    nbins = int(nbins1[0]) # number of range bins
    el1 = elev1 # vertical antenna pointing angle (deg)
    el2 = elev2 # vertical antenna pointing angle (deg)
    el3 = elev3 # vertical antenna pointing angle (deg)
    bw = beamwidth1 # half power beam width (deg)
    range_res = range_res1[0] # range resolution (meters)
    
    r = np.arange(nbins) * range_res
    beamradius = wrl.util.half_power_radius(r, bw)
    
    # elevation 0.5
    
    coord = wrl.georef.sweep_centroids(nrays, range_res, nbins, el1)
    coords = wrl.georef.spherical_to_proj(coord[..., 0],
                                          coord[..., 1],
                                          coord[..., 2], sitecoords)
    lon = coords[..., 0]
    lat = coords[..., 1]
    alt1 = coords[..., 2]
    
    polcoords = coords[..., :2]
    print("lon,lat,alt:", coords.shape)
    
    rlimits = (lon.min(), lat.min(), lon.max(), lat.max())
    print("Radar bounding box:\n\t%.2f\n%.2f             %.2f\n\t%.2f" %
          (lat.max(), lon.min(), lon.max(), lat.min()))
    
    # rasterfile = wrl.util.get_wradlib_data_file('IDN_msk_alt.tif')
    # rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    # ds = wrl.io.open_raster(rasterfile)
    # rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-32768.)
    rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-9999.)
    rastervalues = np.where(rastervalues == -9999, np.nan, rastervalues)
    
    # Clip the region inside our bounding box
    ind = wrl.util.find_bbox_indices(rastercoords, rlimits)
    rastercoords = rastercoords[ind[1]:ind[3], ind[0]:ind[2], ...]
    rastervalues = rastervalues[ind[1]:ind[3], ind[0]:ind[2]]
    
    # Map rastervalues to polar grid points
    polarvalues1 = wrl.ipol.cart_to_irregular_spline(rastercoords, rastervalues,
                                                 polcoords, order=3,
                                                 prefilter=False)
    
    PBB1 = wrl.qual.beam_block_frac(polarvalues1, alt1, beamradius)
    PBB1 = np.ma.masked_invalid(PBB1)
    print(PBB1.shape)
    
    CBB1 = wrl.qual.cum_beam_block_frac(PBB1)
    print(CBB1.shape)
    
    # elevation 1.0
    
    coord = wrl.georef.sweep_centroids(nrays, range_res, nbins, el2)
    coords = wrl.georef.spherical_to_proj(coord[..., 0],
                                          coord[..., 1],
                                          coord[..., 2], sitecoords)
    lon = coords[..., 0]
    lat = coords[..., 1]
    alt2 = coords[..., 2]
    
    polcoords = coords[..., :2]
    print("lon,lat,alt:", coords.shape)
    
    rlimits = (lon.min(), lat.min(), lon.max(), lat.max())
    print("Radar bounding box:\n\t%.2f\n%.2f             %.2f\n\t%.2f" %
          (lat.max(), lon.min(), lon.max(), lat.min()))
    
    # rasterfile = wrl.util.get_wradlib_data_file('IDN_msk_alt.tif')
    # rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    # ds = wrl.io.open_raster(rasterfile)
    # rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-32768.)
    rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-9999.)
    rastervalues = np.where(rastervalues == -9999, np.nan, rastervalues)
    
    # Clip the region inside our bounding box
    ind = wrl.util.find_bbox_indices(rastercoords, rlimits)
    rastercoords = rastercoords[ind[1]:ind[3], ind[0]:ind[2], ...]
    rastervalues = rastervalues[ind[1]:ind[3], ind[0]:ind[2]]
    
    # Map rastervalues to polar grid points
    polarvalues2 = wrl.ipol.cart_to_irregular_spline(rastercoords, rastervalues,
                                                 polcoords, order=3,
                                                 prefilter=False)
    
    PBB2 = wrl.qual.beam_block_frac(polarvalues2, alt2, beamradius)
    PBB2 = np.ma.masked_invalid(PBB2)
    print(PBB2.shape)
    
    CBB2 = wrl.qual.cum_beam_block_frac(PBB2)
    print(CBB2.shape)
    
    # elevation 1.5
    
    coord = wrl.georef.sweep_centroids(nrays, range_res, nbins, el3)
    coords = wrl.georef.spherical_to_proj(coord[..., 0],
                                          coord[..., 1],
                                          coord[..., 2], sitecoords)
    lon = coords[..., 0]
    lat = coords[..., 1]
    alt3 = coords[..., 2]
    
    polcoords = coords[..., :2]
    print("lon,lat,alt:", coords.shape)
    
    rlimits = (lon.min(), lat.min(), lon.max(), lat.max())
    print("Radar bounding box:\n\t%.2f\n%.2f             %.2f\n\t%.2f" %
          (lat.max(), lon.min(), lon.max(), lat.min()))
    
    # rasterfile = wrl.util.get_wradlib_data_file('IDN_msk_alt.tif')
    # rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    # ds = wrl.io.open_raster(rasterfile)
    # rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-32768.)
    rastervalues, rastercoords, proj = wrl.georef.extract_raster_dataset(ds, nodata=-9999.)
    rastervalues = np.where(rastervalues == -9999, np.nan, rastervalues)
    
    # Clip the region inside our bounding box
    ind = wrl.util.find_bbox_indices(rastercoords, rlimits)
    rastercoords = rastercoords[ind[1]:ind[3], ind[0]:ind[2], ...]
    rastervalues = rastervalues[ind[1]:ind[3], ind[0]:ind[2]]
    
    # Map rastervalues to polar grid points
    polarvalues3 = wrl.ipol.cart_to_irregular_spline(rastercoords, rastervalues,
                                                 polcoords, order=3,
                                                 prefilter=False)
    
    PBB3 = wrl.qual.beam_block_frac(polarvalues3, alt3, beamradius)
    PBB3 = np.ma.masked_invalid(PBB3)
    print(PBB3.shape)
    
    CBB3 = wrl.qual.cum_beam_block_frac(PBB3)
    print(CBB3.shape)
    
    
    fig = pl.figure(figsize=(10, 8))
    
    # create subplots
    ax1 = pl.subplot2grid((3, 2), (0, 0))
    ax2 = pl.subplot2grid((3, 2), (0, 1))
    ax3 = pl.subplot2grid((3, 2), (1, 0))
    ax4 = pl.subplot2grid((3, 2), (1, 1))
    ax5 = pl.subplot2grid((3, 2), (2, 0), colspan=2, rowspan=1)
    
    # azimuth angle
    angle = azi_angle
    
    # Plot terrain (on ax1)
    ax1, dem = wrl.vis.plot_ppi(polarvalues1,
                                ax=ax1, r=r,
                                az=coord[:,0,1],
                                cmap=mpl.cm.terrain, vmin=0.)
    ax1.plot([0,np.sin(np.radians(angle))*r.max()],
             [0,np.cos(np.radians(angle))*r.max()],"r-")
    ax1.plot(sitecoords[0], sitecoords[1], 'ro')
    annotate_map(ax1, dem, 'Terrain within {:.0f} km range'.format(np.max(r / 1000.) + range_res/1000.) + ' (' + sitename1+')')
    
    # Plot CBB (on ax2)
    ax2, cbb = wrl.vis.plot_ppi(CBB1, ax=ax2, r=r,
                                az=coord[:,0,1],
                                cmap=mpl.cm.PuRd, vmin=0, vmax=1)
    annotate_map(ax2, cbb, 'Beam-Blockage Fraction ' + np.str(el1) +' deg (' + sitename1+')')
    
    # Plot CBB (on ax3)
    ax3, cbb = wrl.vis.plot_ppi(CBB2, ax=ax3, r=r,
                                az=coord[:,0,1],
                                cmap=mpl.cm.PuRd, vmin=0, vmax=1)
    annotate_map(ax3, cbb, 'Beam-Blockage Fraction ' + np.str(el2) +' deg (' + sitename1+')')
    
    # Plot CBB (on ax4)
    ax4, cbb = wrl.vis.plot_ppi(CBB3, ax=ax4, r=r,
                                az=coord[:,0,1],
                                cmap=mpl.cm.PuRd, vmin=0, vmax=1)
    annotate_map(ax4, cbb, 'Beam-Blockage Fraction ' + np.str(el3) +' deg (' + sitename1+')')
    
    
    # Plot single ray terrain profile on ax3
    bc1, = ax5.plot(r / 1000., alt1[angle, :], '-r',
                   linewidth=2, label='Beam Center (' + str(el1) + ' deg)')
    bc2, = ax5.plot(r / 1000., alt2[angle, :], '-b',
                   linewidth=2, label='Beam Center (' + str(el2) + ' deg)')
    bc3, = ax5.plot(r / 1000., alt3[angle, :], '-g',
                   linewidth=2, label='Beam Center (' + str(el3) + ' deg)')
    # b3db, = ax5.plot(r / 1000., (alt2[angle, :] + beamradius), ':b',
    #                  linewidth=1.5, label='3 dB Beam width')
    # ax5.plot(r / 1000., (alt2[angle, :] - beamradius), ':b')
    ax5.fill_between(r / 1000., 0.,
                     polarvalues1[angle, :],
                     color='0.75')
    ax5.set_xlim(0., np.max(r / 1000.) + range_res/1000.)
    ax5.set_ylim(0., 5000)
    ax5.set_xlabel('Range (km)')
    ax5.set_ylabel('Altitude (m)')
    ax5.grid()
    
    axb = ax5.twinx()
    bbf, = axb.plot(r / 1000., CBB2[angle, :], '-k',
                    label='BBF')
    axb.set_ylabel('Beam-blockage fraction')
    axb.set_ylim(0., 1.)
    axb.set_xlim(0., np.max(r / 1000.) + range_res/1000.)
    
    legend = ax5.legend((bc1, bc2, bc3, bbf),
                        ('Beam Center ('+ str(el1) + ' deg)', 'Beam Center ('+ str(el2) + ' deg)', 'Beam Center ('+ str(el3) + ' deg)', 'BBF'),
                        loc='upper left', fontsize=10)
    # legend = ax5.legend((bc1, bc2, bc3, b3db, bbf),
    #                     ('Beam Center ('+ str(el1) + ' deg)', 'Beam Center ('+ str(el2) + ' deg)', 'Beam Center ('+ str(el3) + ' deg)', '3 dB Beam width', 'BBF'),
    #                     loc='upper left', fontsize=10)
    
    pl.annotate('DSP',xy=(0.995, -0.30), xycoords='axes fraction', fontsize=5,horizontalalignment='right', verticalalignment='bottom')
    
    # Save Figure
    main_directory = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/'
    # pl.savefig(main_directory+'BBF_'+sitename1+'.pdf',format='pdf',dpi=100,bbox_inches='tight')
    # pl.clf()
    pl.savefig(main_directory+'BBF_'+sitename1+'.png',dpi=300,bbox_inches='tight')
    
filenamelist1 = []
azi_anglelist1 = []

filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\ACE2020080700000600dBZ.vol')
azi_anglelist1.append( [165, 18]) # azimuth angle, tower height
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\BKS2020080700000700dBZ.vol')
azi_anglelist1.append( [30, 20])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\BMA2020080700100500dBZ.vol')
azi_anglelist1.append( [250, 17])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\DJB2020080700193500dBZ.vol')
azi_anglelist1.append( [220, 23])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\LOP2020080700000600dBZ.vol')
azi_anglelist1.append( [0, 20])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\MAU2020080723500200dBZ.vol')
azi_anglelist1.append( [310, 17])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\MND2020080700000200dBZ.vol')
azi_anglelist1.append( [0, 12])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\MRS2019041800001300dBZ.vol')
azi_anglelist1.append( [80, 19])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\PDG2020080700062300dBZ.vol')
azi_anglelist1.append( [30, 20])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\PKN2020080700000400dBZ.vol')
azi_anglelist1.append( [0, 20])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\PKY2020080600000500dBZ.vol')
azi_anglelist1.append( [100, 24])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\PLB2020080700000500dBZ.vol')
azi_anglelist1.append( [220, 20])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\SBY2020080700000600dBZ.vol')
azi_anglelist1.append( [170, 23])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\SIN2020080702035500dBZ.vol')
azi_anglelist1.append( [205, 21])
filenamelist1.append('I:\Radar_Integrasi\\2020\GEMA\TTE2020080700000600dBZ.vol')
azi_anglelist1.append( [125, 18])

h = 0
for filename in filenamelist1: 
    azi_angle = azi_anglelist1[h][0]
    tower_height = azi_anglelist1[h][1]
    h = h + 1
    rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    ds = wrl.io.open_raster(rasterfile)
    
    # raw = wrl.io.read_Rainbow(filename)
    sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename = process_radar_data_gematronik(filename)
    sitecoords = tuple(map(sum, zip(sitecoords, (0,0,tower_height))))    
    plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)



filenamelist2 = []
azi_anglelist2 = []

filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\JAK-20200501-140214-PPIVol.nc')
azi_anglelist2.append( [170, 20])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\1004KEN-20200501-184002-PPIVol.nc')
azi_anglelist2.append( [300, 16])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\1013SOR-20200501-181001-PPIVol.nc')
azi_anglelist2.append( [130, 20])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\1200TIM-20200212-185202-PPIVol.nc')
azi_anglelist2.append( [0, 20])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\9257GOR-20200501-170429-PPIVol.nc')
azi_anglelist2.append( [170, 20])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\9262BAN-20200501-170002-PPIVol.nc')
azi_anglelist2.append( [130, 23])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\\9501TAR-20200308-025002-PPIVol.nc')
azi_anglelist2.append( [260, 16])
# filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\BAL-20200501-141006-PPIVol.nc')
# azi_anglelist2.append( [170, 25])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\BIK-20200501-180004-PPIVol.nc')
azi_anglelist2.append( [260, 24])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\DEN-20200501-175005-PPIVol.nc')
azi_anglelist2.append( [0, 25])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\LPG-20200501-174005-PPIVol.nc')
azi_anglelist2.append( [270, 27])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\MED-20200501-182006-PPIVol.nc')
azi_anglelist2.append( [260, 23])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\PEK-20200501-184542-PPIVol.nc')
azi_anglelist2.append( [215, 33])
filenamelist2.append('I:\Radar_Integrasi\\2020\EEC\PON-20200501-180436-PPIVol.nc')
azi_anglelist2.append( [10, 26])

h = 0
for filename in filenamelist2: 
    azi_angle = azi_anglelist2[h][0]
    tower_height = azi_anglelist2[h][1]
    h = h + 1
    rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    ds = wrl.io.open_raster(rasterfile)
    
    # raw = wrl.io.read_Rainbow(filename)
    sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename = process_radar_data_eec(filename)
    sitecoords = tuple(map(sum, zip(sitecoords, (0,0,tower_height))))
    plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)


filenamelist3 = []
azi_anglelist3 = []

filenamelist3.append('I:\Radar_Integrasi\\2020\BARON\GNS20181121035004.vol')
azi_anglelist3.append( [170, 16])
filenamelist3.append('I:\Radar_Integrasi\\2020\BARON\KPG20181211041000.vol')
azi_anglelist3.append( [90, 15])
filenamelist3.append('I:\Radar_Integrasi\\2020\BARON\PKP20200202142145.vol')
azi_anglelist3.append( [200, 18])
filenamelist3.append('I:\Radar_Integrasi\\2020\BARON\SMG20200303002000.vol')
azi_anglelist3.append( [135, 15])
filenamelist3.append('I:\Radar_Integrasi\\2020\BARON\YOG20200303015404.vol')
azi_anglelist3.append( [270, 16])

h = 0
for filename in filenamelist3: 
    azi_angle = azi_anglelist3[h][0]
    tower_height = azi_anglelist3[h][1]
    h = h + 1
    rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
    
    ds = wrl.io.open_raster(rasterfile)
    
    # raw = wrl.io.read_Rainbow(filename)
    sitecoords, elevation, nrays, nbins, range_res, beamwidth, sitename = process_radar_data_baron(filename)
    sitecoords = tuple(map(sum, zip(sitecoords, (0,0,tower_height))))
    plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, os.path.basename(filename), ds, 0.5, 1.0, 1.5)


rasterfile = 'G:/DATA/BMKG/2020_BMKG_GAW/08/08_Radar_Coverage/IDN_msk_alt.tif'
ds = wrl.io.open_raster(rasterfile)

sitename = 'PT. Freeport Indonesia'
sitecoords = (137.118236, -4.046066, 4000 + 15)
nrays = [360] # number of rays
nbins = [800] # number of range bins # C-Band ~ 200 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 180
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, -0.5, 0.0, 0.5)

sitename = 'Jatiluhur'
sitecoords = (107.417445, -6.530567, 114 + 20)
nrays = [360] # number of rays
nbins = [200] # number of range bins # X-Band ~ 50 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 180
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'Lembang'
sitecoords = (107.617746, -6.826462, 1253 + 15)
nrays = [360] # number of rays
nbins = [800] # number of range bins # C-Band ~ 200 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 270
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'Sragen'
sitecoords = (111.023394, -7.426977, 93 + 15)
nrays = [360] # number of rays
nbins = [200] # number of range bins # X-Band ~ 50 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 135
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'Nabire'
sitecoords = (135.498037, -3.36594, 12 + 15)
nrays = [360] # number of rays
nbins = [800] # number of range bins # C-Band ~ 200 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 90
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'Pulau Obi(Malut)'
sitecoords = (127.637895, -1.340005, 9 + 15)
nrays = [360] # number of rays
nbins = [800] # number of range bins # C-Band ~ 200 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 0
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'Morotai'
sitecoords = (128.322494, 2.042035, 17 + 15)
nrays = [360] # number of rays
nbins = [200] # number of range bins # X-Band ~ 50 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 30
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)

sitename = 'PT Gudang Garam'
sitecoords = (111.950193, -7.774645, 200 + 15)
nrays = [360] # number of rays
nbins = [200] # number of range bins # X-Band ~ 50 km
# el = 1.0 #elevation1[1] # vertical antenna pointing angle (deg)
beamwidth = 1 # half power beam width (deg)
range_res = [250] # range resolution (meters)
azi_angle = 270
plotBBFanalysis(sitecoords, nrays, nbins, beamwidth, range_res, sitename, ds, 0.5, 1.0, 1.5)