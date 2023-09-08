#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as image
import mpl_toolkits
mpl_toolkits.__path__.append('/gpfs/dell2/emc/modeling/noscrub/gwv/py/lib/python/basemap-1.2.1-py3.6-linux-x86_64.egg/mpl_toolkits/')
from mpl_toolkits.basemap import Basemap, cm
#from matplotlib import GridSpec, rcParams, colors
from matplotlib import colors as c

import netCDF4
import numpy as np
import pygrib, datetime, time, os, sys, subprocess
import ncepy, scipy
from ncepgrib2 import Grib2Encode, Grib2Decode
import dawsonpy
import itertools

def cmap_t2m():
 # Create colormap for 2-m temperature
 # Modified version of the ncl_t2m colormap from Jacob's ncepy code
    r=np.array([255,128,0,  70, 51, 0,  255,0, 0,  51, 255,255,255,255,255,171,128,128,36,162,255])
    g=np.array([0,  0,  0,  70, 102,162,255,92,128,185,255,214,153,102,0,  0,  0,  68, 36,162,255])
    b=np.array([255,128,128,255,255,255,255,0, 0,  102,0,  112,0,  0,  0,  56, 0,  68, 36,162,255])
    xsize=np.arange(np.size(r))
    r = r/255.
    g = g/255.
    b = b/255.
    red = []
    green = []
    blue = []
    for i in range(len(xsize)):
        xNorm=np.float(i)/(np.float(np.size(r))-1.0)
        red.append([xNorm,r[i],r[i]])
        green.append([xNorm,g[i],g[i]])
        blue.append([xNorm,b[i],b[i]])
    colorDict = {"red":red, "green":green, "blue":blue}
    cmap_t2m_coltbl = matplotlib.colors.LinearSegmentedColormap('CMAP_T2M_COLTBL',colorDict)
    return cmap_t2m_coltbl


# Get machine and head directory
machine, hostname = dawsonpy.get_machine()

case = 'derecho'
machine = 'WCOSS_DELL_P3'
if machine == 'WCOSS':
    DIR = '/gpfs/'+hostname[0]+'p2/ptmp/'+os.environ['USER']+'/rap_hrrr_retros'
    DATA_DIR = os.path.join('/gpfs/dell2/emc/verification/noscrub/Geoffrey.Manikin/',case)
    GRAPHX_DIR = os.path.join('/gpfs/dell1/ptmp/Geoffrey.Manikin/'+case+'graphx')
elif machine == 'WCOSS_DELL_P3':
    DATA_DIR = os.path.join('/gpfs/dell2/emc/verification/noscrub/Geoffrey.Manikin/',case)
    GRAPHX_DIR = os.path.join('/gpfs/dell1/ptmp/Geoffrey.Manikin/'+case+'graphx')


# Create graphx directory
if os.path.exists(DATA_DIR):
   if not os.path.exists(GRAPHX_DIR):
      os.makedirs(GRAPHX_DIR)
else:
#  raise NameError, 'data for '+case+' case not found'
   print('data for '+case+' case not found. Graphics will be saved in current directory.')
   DIR = os.getcwd()
   GRAPHX_DIR = os.path.join(DIR,'graphx')
   if not os.path.exists(GRAPHX_DIR):
      os.makedirs(GRAPHX_DIR)


# Determine initial date/time
try:
   cycle = str(sys.argv[1])
except IndexError:
   cycle = None

if cycle is None:
   cycle = input('Enter initial time (YYYYMMDDHH): ')

YYYY = int(cycle[0:4])
MM   = int(cycle[4:6])
DD   = int(cycle[6:8])
HH   = int(cycle[8:10])
print(YYYY, MM, DD, HH)

date_str = datetime.datetime(YYYY,MM,DD,HH,0,0)


#Build forecast hour list
if cycle == '2020080912':
   fhrs = np.arange(19,37,1)
elif cycle == '2020080918':
   fhrs = np.arange(13,34,1)

elif cycle == '2020081000':
   fhrs = np.arange(7,28,1)
elif cycle == '2020081001':
   fhrs = np.arange(6,27,1)
elif cycle == '2020081002':
   fhrs = np.arange(5,26,1)
elif cycle == '2020081003':
   fhrs = np.arange(4,19,1)
elif cycle == '2020081004':
   fhrs = np.arange(3,19,1)
elif cycle == '2020081005':
   fhrs = np.arange(2,19,1)
elif cycle == '2020081006':
   fhrs = np.arange(1,22,1)

elif cycle == '2020081007':
   fhrs = np.arange(1,19,1)
elif cycle == '2020081008':
   fhrs = np.arange(1,19,1)
elif cycle == '2020081009':
   fhrs = np.arange(1,19,1)
elif cycle == '2020081010':
   fhrs = np.arange(1,18,1)
elif cycle == '2020081011':
   fhrs = np.arange(1,17,1)
elif cycle == '2020081012':
   fhrs = np.arange(1,16,1)

elif cycle == '2020081013':
   fhrs = np.arange(1,15,1)
elif cycle == '2020081014':
   fhrs = np.arange(1,14,1)
elif cycle == '2020081015':
   fhrs = np.arange(1,13,1)
elif cycle == '2020081016':
   fhrs = np.arange(1,12,1)
elif cycle == '2020081017':
   fhrs = np.arange(1,11,1)

valid_list = [date_str + datetime.timedelta(hours=int(x)) for x in fhrs]


large_domains = ['CONUS','conus','Australia']
regional_domains = ['eCONUS']
subregional_domains = ['SPL','MIDSOUTH','Barry']
state_domains = ['OK','Louisiana','TN']


#Specify plots to make
fields = ['t2_10mwind','refc_uh','sbcape','mlcape']
field = fields[1]
field = 'max10mwind'
#field = 'refc'
domain = 'derecho'

maxwind_list = []

for i in range(1):
   print(valid_list[0],valid_list[-1])
   valid = valid_list[-1].strftime("%Y%m%d%H")

   fig = plt.figure(figsize=(8,8))
   ax1 = fig.add_subplot(111)
   axes = [ax1]

   k = 0
   for ax in axes:

      print('plotting '+str.lower(field)+' from F'+str(fhrs[0]-1).zfill(2)+' to F'+str(fhrs[-1]).zfill(2)+' on '+domain+' domain')

      if domain == 'conus':
          print('specify new corner vals for CONUS')
          llcrnrlon=-121.0
          llcrnrlat=21.0
          urcrnrlon=-62.6
          urcrnrlat=49.0
      elif domain == 'nw':
          llcrnrlat=36.
          llcrnrlon=-126.
          urcrnrlat=53.
          urcrnrlon=-108.
      elif domain == 'nc':
          llcrnrlat=36.
          llcrnrlon=-112.
          urcrnrlat=52.5
          urcrnrlon=-82.
      elif domain == 'ne':
          llcrnrlat=39.
          llcrnrlon=-89.
          urcrnrlat=49.
          urcrnrlon=-63.
      elif domain == 'sw':
          llcrnrlat=22.
          llcrnrlon=-122.
          urcrnrlat=41.
          urcrnrlon=-106.
      elif domain == 'sc':
          llcrnrlat=24.
          llcrnrlon=-109.
          urcrnrlat=41.
          urcrnrlon=-85.
      elif domain == 'se':
          llcrnrlat=24.
          llcrnrlon=-91.
          urcrnrlat=38.
          urcrnrlon=-68.
      elif domain == 'midatl':
          llcrnrlat=34.0
          llcrnrlon=-85.5
          urcrnrlat=41.0
          urcrnrlon=-71.25
      elif domain == 'glakes':
          llcrnrlat=40.5
          llcrnrlon=-93.5
          urcrnrlat=48.5
          urcrnrlon=-74.2
      elif domain == 'ca':
          llcrnrlat=32.
          llcrnrlon=-124.
          urcrnrlat=42.
          urcrnrlon=-117.
      elif domain == 'derecho':
          llcrnrlat=36.
          llcrnrlon=-102.
          urcrnrlat=46.
          urcrnrlon=-82.

      m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                  resolution='i',projection='lcc',\
                  lat_1=32.,lat_2=46.,lon_0=-101.,area_thresh=1000.,ax=ax)


      m.drawcoastlines()
      m.drawstates(linewidth=0.75)
      m.drawcountries()

      if domain in large_domains:
         latlongrid = 10.
         barb_length = 4.5
         parallels = np.arange(-90.,91.,latlongrid)
         m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
         meridians = np.arange(0.,360.,latlongrid)
         m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
      else:
         m.drawcounties(zorder=20)
         latlongrid = 5.
         barb_length = 5.5


      # set this as default, and will update if model doesn't exist for given cycle
      model_cycle = 'exists'

      if k == 0:
         model_str = 'NAM NEST'

         # Loop over forecast hours and append max wind data to a list
         for j in range(len(valid_list)):
             print('reading F'+str(fhrs[j]).zfill(2))
             fil = DATA_DIR+'/nam.'+cycle[0:8]+'/nam.t'+ \
                   cycle[8:10]+'z.conusnest.hiresf'+str(fhrs[j]).zfill(2)+'.tm00.grib2'
             grbs = pygrib.open(fil)

             u_var = grbs.select(parameterName='222',typeOfLevel='heightAboveGround')[0].values*1.94384
             v_var = grbs.select(parameterName='223',typeOfLevel='heightAboveGround')[0].values*1.94384
             max_wind = np.sqrt(u_var**2+v_var**2)

             maxwind_list.append(max_wind)

             # Close grib file for each fhr except the last
             # Leave last one open to read lat/losn below
             if fhrs[j] != fhrs[-1]:
                 grbs.close()

         # After loop, take the max value at each gridpoint
         maxwind_accum = np.amax(np.array(maxwind_list),axis=0)
         

      elif k == 1:
         model_str = 'HiResWARW'
         try:
             fil = DATA_DIR+'/href.'+cycle[0:8]+'/'+str.lower(model_str[0:6])+'.t'+ \
                   cycle[8:10]+'z.'+str.lower(model_str[6:])+'_5km.f'+str(fhrs[j]).zfill(2)+'.conus.grib2'
             grbs = pygrib.open(fil)
         except:
             model_cycle = None

      elif k == 2:
         model_str = 'HiResWARW2'
         try:
             fil = DATA_DIR+'/href.'+cycle[0:8]+'/'+str.lower(model_str[0:6])+'.t'+ \
                   cycle[8:10]+'z.'+str.lower(model_str[6:9])+'_5km.f'+str(fhrs[j]).zfill(2)+'.conusmem2.grib2'
             grbs = pygrib.open(fil)
         except:
             model_cycle = None

      elif k == 3:
         model_str = 'NAM3'

         # Loop over forecast hours and append max wind data to a list
         for j in range(len(valid_list)):
             print('reading F'+str(fhrs[j]).zfill(2))
             fil = DATA_DIR+'/nam.'+cycle[0:8]+'/'+str.lower(model_str[0:3])+'.t'+ \
                   cycle[8:10]+'z.conusnest.hiresf'+str(fhrs[j]).zfill(2)+'.tm00.grib2'
             grbs = pygrib.open(fil)

             u_var = grbs.select(parameterName='222',typeOfLevel='heightAboveGround')[0].values*1.94384
             v_var = grbs.select(parameterName='223',typeOfLevel='heightAboveGround')[0].values*1.94384
             max_wind = np.sqrt(u_var**2+v_var**2)

             maxwind_list.append(max_wind)

             # Close grib file for each fhr except the last
             # Leave last one open to read lat/losn below
             if fhrs[j] != fhrs[-1]:
                 grbs.close()

         # After loop, take the max value at each gridpoint
         maxwind_accum = np.amax(np.array(maxwind_list),axis=0)
         

      elif k == 4:
         model_str = 'HiResWNMMB'
         try:
             fil = DATA_DIR+'/href.'+cycle[0:8]+'/'+str.lower(model_str[0:6])+'.t'+ \
                   cycle[8:10]+'z.'+str.lower(model_str[6:])+'_5km.f'+str(fhrs[j]).zfill(2)+'.conus.grib2'
             grbs = pygrib.open(fil)
         except:
             model_cycle = None


      # Obs panel - MRMS for REFC
      elif k == 5 and field[0:4] == 'refc':
         model_str = 'MRMS Observations'
         fil = '/gpfs/dell2/ptmp/Logan.Dawson/com/MRMS/derecho/'+ \
               'MergedReflectivityQCComposite_00.50_'+valid[0:8]+'-'+valid[8:10]+'0000.nc'
         nc = netCDF4.Dataset(fil,'r')
         refc = nc.variables['MergedReflectivityQCComposite'][:]
         grid_lats = nc.variables['lat'][:]
         grid_lons = nc.variables['lon'][:]
         nc.close()


      # Obs panel - URMA for 2-m T
      elif k == 5 and field[0:2] == 't2':
         model_str = 'URMA Analysis'
         anl_str = 'URMA'
         try:
             fil = '/gpfs/dell2/nco/ops/com/'+str.lower(anl_str)+'/prod/'+str.lower(anl_str)+'2p5.'+valid[0:8]+'/'+ \
                   str.lower(anl_str)+'2p5.t'+valid[8:10]+'z.2dvaranl_ndfd.grb2_wexp'
             grbs = pygrib.open(fil)
         except:
             fil = DATA_DIR+'/'+str.lower(anl_str)+'2p5.'+valid[0:8]+'.t'+ \
                   valid[8:10]+'z.2dvaranl_ndfd.grb2'
             grbs = pygrib.open(fil)

      # Obs panel - RAP for CAPE
      elif k == 5 and field[2:6] == 'cape':
         model_str = 'RAP Analysis'
         anl_str = 'RAP'
         try:
             fil = '/gpfs/dell2/nco/ops/com/'+str.lower(anl_str)+'/prod/'+str.lower(anl_str)+'2p5.'+valid[0:8]+'/'+ \
                   str.lower(anl_str)+'2p5.t'+valid[8:10]+'z.2dvaranl_ndfd.grb2_wexp'
             grbs = pygrib.open(fil)
         except:
             fil = DATA_DIR+'/'+str.lower(anl_str)+'.'+valid[0:8]+'.t'+ \
                   valid[8:10]+'z.awp130pgrbf00.nc'
             nc = netCDF4.Dataset(fil,'r')
             sbcape = nc.variables['CAPE_L0'][:]
             sbcin  = nc.variables['CIN_L0'][:]
             mlcape = nc.variables['CAPE_P90-0'][:]
             mlcin  = nc.variables['CIN_P90-0'][:]
             ylats  = nc.variables['lat'][:]
             xlons  = nc.variables['lon'][:]
             nc.close()


      # If it's a cycle when HRWs are missing
      # Read a HRRRX file, and we'll just zero it out below
      if model_cycle is None:
         fil = DATA_DIR+'/nam.'+cycle[0:8]+'/nam.t'+ \
                   cycle[8:10]+'z.conusnest.hiresf'+str(fhrs[j]).zfill(2)+'.tm00.grib2'
         print('using this file to zero out: '+fil)
         grbs = pygrib.open(fil)


      if k == 5 and field[0:4] == 'refc':
         print('skipping grib lat/lons bc netcdf data')
      elif k == 5 and field[2:6] == 'cape':
         print('skipping grib lat/lons bc netcdf data')
      elif k == 5 and field == 'max10mwind':
         print('skipping grib lat/lons bc not reading data')
      else:
         grb = grbs.message(1)
         lats, lons = grb.latlons()         



      if str.upper(model_str[0:6]) == 'HIRESW':
         HLwindow = 250
         if domain in large_domains:
            skip = 30
         elif domain in subregional_domains:
            skip = 15
         elif domain in state_domains:
            skip = 10

      elif str.upper(model_str) == 'NAM3' or str.upper(model_str) == 'NAM NEST':
         HLwindow = 300
         if domain in large_domains:
            skip = 50
         elif domain in subregional_domains:
            skip = 25
         elif domain in state_domains:
            skip = 20

      elif str.upper(model_str[0:4]) == 'RTMA' or str.upper(model_str[0:4]) == 'URMA':
         HLwindow = 400
         if domain in large_domains:
            skip = 75
         elif domain in subregional_domains:
            skip = 30
         elif domain in state_domains:
            skip = 20



      # Temperature
      if field[0:2] == 't2':
         t2  = grbs.select(stepType='instant',name='2 metre temperature')[0].values*9/5 - 459.67
         u10 = grbs.select(stepType='instant',name='10 metre U wind component')[0].values*1.94384
         v10 = grbs.select(stepType='instant',name='10 metre V wind component')[0].values*1.94384

         clevs = np.arange(-16,132,4)
         tlevs = [str(clev) for clev in clevs]

         colormap = cmap_t2m()
         norm = matplotlib.colors.BoundaryNorm(clevs, colormap.N)

         fill_var = t2

         # Zero out fill_var if model cycle doesn't exist
         if model_cycle is None:
             fill_var = fill_var*0.
             print(model_str+' doest exist for this '+cycle+' cycle. Zeroing out data.')

         fill = m.contourf(lons,lats,fill_var,clevs,latlon=True,cmap=colormap,norm=norm,extend='both')
         field_title = '2-m Temperature'
         titlestr = date_str.strftime('%HZ %m/%d ')+field_title+' Forecasts and Analysis\n'+valid_list[j].strftime('Valid at %HZ %d %b %Y')+' (F'+str(fhrs[j]).zfill(2)+')'

         try:
            if str.lower(field[3:10]) == '10mwind':

               u_var = u10
               v_var = v10

               m.barbs(lons[::skip,::skip],lats[::skip,::skip],u_var[::skip,::skip],v_var[::skip,::skip],latlon=True,length=barb_length,sizes={'spacing':0.2},pivot='middle')

            #  field_title = '2-m Temperature and 10-m Wind'

         except IndexError:
            sys.exc_clear()


      # CAPE
      elif field[2:6] == 'cape':
         print(fil)
         if k != 5:
            if field[0:2] == 'sb':
               sbcape = grbs.select(name='Convective available potential energy',typeOfLevel='surface')[0].values
               sbcin  = grbs.select(name='Convective inhibition',typeOfLevel='surface')[0].values

            elif field[0:2] == 'ml':
               mlcape = grbs.select(name='Convective available potential energy',typeOfLevel='pressureFromGroundLayer',topLevel=9000)[0].values
               mlcin  = grbs.select(name='Convective inhibition',typeOfLevel='pressureFromGroundLayer',topLevel=9000)[0].values

         clevs = [100,250,500,750,1000,1250,1500,1750,2000,2500,3000,3500,4000,5000]
         tlevs = [int(clev) for clev in clevs]
         colorlist = ['lightgray','silver','darkgray','gray', \
                      'lightblue','skyblue','cornflowerblue','steelblue', \
                      'chartreuse','limegreen','yellow','gold','darkorange','red']

         if str.lower(field[0:2]) == 'sb':
            fill_var = sbcape
         elif str.lower(field[0:2]) == 'ml':
            fill_var = mlcape
         elif str.lower(field[0:2]) == 'mu':
            fill_var = mucape

         fill_var = scipy.ndimage.gaussian_filter(fill_var,1)
         if k != 5:
            fill_var = scipy.ndimage.gaussian_filter(fill_var,1)

         # Zero out fill_var if model cycle doesn't exist
         if model_cycle is None:
             fill_var = fill_var*0.
             print(model_str+' doest exist for this '+cycle+' cycle. Zeroing out data.')

         if k != 5:
            fill = m.contourf(lons,lats,fill_var,clevs,latlon=True,colors=colorlist,extend='max')
         else:
            fill = m.contourf(xlons,ylats,fill_var,clevs,latlon=True,colors=colorlist,extend='max')

         field_title = str.upper(field[0:6])
         titlestr = date_str.strftime('%HZ %m/%d ')+field_title+' Forecasts and Analysis\n'+valid_list[j].strftime('Valid at %HZ %d %b %Y')+' (F'+str(fhrs[j]).zfill(2)+')'

         try:
            if str.lower(field[0:2]) == 'sb':
               contour_var = sbcin
            elif str.lower(field[0:2]) == 'ml':
               contour_var = mlcin

            contour_var = scipy.ndimage.gaussian_filter(contour_var,2)
            if k != 5:
               contour_var = scipy.ndimage.gaussian_filter(contour_var,2)


            cint = [-100.,-25.]
            if k != 5:
               contours = m.contour(lons,lats,contour_var,cint,colors='k',linewidths=0.5,linestyles='solid',latlon=True)
               hatch = m.contourf(lons,lats,contour_var,cint,extend='min',colors='none',hatches=['\/\/','..'],latlon=True)
            else:
               contours = m.contour(xlons,ylats,contour_var,cint,colors='k',linewidths=0.5,linestyles='solid',latlon=True)
               hatch = m.contourf(xlons,ylats,contour_var,cint,extend='min',colors='none',hatches=['\/\/','..'],latlon=True)

#           plt.clabel(contours,cint,colors='k',inline=1,fmt='%.0f',fontsize=9)

         except IndexError:
            sys.exc_clear()


      # Simulated reflectivity
      elif field[0:3] == 'ref':
         if k != 5:
             refc = grbs.select(name='Maximum/Composite radar reflectivity')[0].values
      
         clevs = np.linspace(5,70,14)
         tlevs = [int(clev) for clev in clevs]
         colorlist = ['turquoise','dodgerblue','mediumblue','lime','limegreen','green', \
                      'yellow','gold','darkorange','red','firebrick','darkred','fuchsia','darkmagenta']

         if field[0:4] == 'refc':
            fill_var = refc
            field_title = 'Composite Reflectivity'
         elif field[0:5] == 'refd1':
            fill_var = refd1
         elif field[0:5] == 'refd4':
            fill_var = refd4

         # Zero out fill_var if model cycle doesn't exist
         if model_cycle is None:
             fill_var = fill_var*0.
             print(model_str+' doest exist for this '+cycle+' cycle. Zeroing out data.')

         if k != 5:
            fill = m.contourf(lons,lats,fill_var,clevs,latlon=True,colors=colorlist,extend='max')
         else:
         #  continue
         #  lons, lats = np.meshgrid(grid_lons,grid_lats)
         #  xlons, ylats = m(lons,lats)
         #  fill = m.contourf(xlons,ylats,fill_var,clevs,colors=colorlist,extend='max')
            fill = m.contourf(grid_lons,grid_lats,fill_var,clevs,colors=colorlist,latlon=True,extend='max')

         titlestr = date_str.strftime('%HZ %m/%d ')+field_title+' Forecasts and Observations\n'+valid_list[j].strftime('Valid at %HZ %d %b %Y')+' (F'+str(fhrs[j]).zfill(2)+')'

      # Hourly-max 10-m wind
      elif field == 'max10mwind':
         if k != 5:
             fill_var = maxwind_accum
             accum_length = fhrs[-1] - fhrs[0] + 1
             field_title = str(accum_length)+'-h Max 10-m Wind Speed'

         clevs = np.arange(15.,75.,5.)
         tlevs = [int(clev) for clev in clevs]
         cmap = plt.get_cmap(name='rainbow')

         # Zero out fill_var if model cycle doesn't exist
         if model_cycle is None:
             fill_var = fill_var*0.
             print(model_str+' doest exist for this '+cycle+' cycle. Zeroing out data.')

         if k != 5:
             fill = m.contourf(lons,lats,fill_var,clevs,latlon=True,cmap=cmap,extend='max')
             fill.set_clim(15,70)
         else:
             continue

         titlestr = date_str.strftime('%HZ %m/%d ')+field_title+' Forecasts\n'+valid_list[j].strftime('Valid at %HZ %d %b %Y')+' (F'+str(fhrs[j]).zfill(2)+')'



      if str.upper(model_str) == 'NAM3':
          mod_str = str.upper(model_str[0:3])+' Nest'
      elif str.upper(model_str[0:6]) == 'HIRESW':
          mod_str = model_str[0:6]+' '+str.upper(model_str[6:])
      elif str.lower(model_str[4:]) == 'mean':
          mod_str = str.upper(model_str[0:4])+' Mean'
      else:
          mod_str = str.upper(model_str)
          mod_str = model_str

      ax.text(0.5, 0.05, mod_str,  horizontalalignment='center', weight='bold', transform=ax.transAxes, bbox=dict(facecolor='white',alpha=0.85))

      if k == 5 and field[0:4] == 'refc':
         print('skipping grib close bc netcdf data')
      elif k == 5 and field[2:6] == 'cape':
         print('skipping grib close bc netcdf data')
      elif k == 5 and field == 'max10mwind':
         print('skipping grib close bc not reading data')
      else:
         grbs.close()
      k += 1


   cax = fig.add_axes([0.2,0.01,0.6,0.03])
   cbar = plt.colorbar(fill,cax=cax,ticks=clevs,orientation='horizontal')
   cbar.ax.set_xticklabels(tlevs)
   cbar.ax.tick_params(labelsize=10)

   for label in cbar.ax.xaxis.get_ticklabels():
      label.set_visible(False)

   if field[0:2] == 't2':
      cbar.set_label(u'\xb0''F')

      for label in cbar.ax.xaxis.get_ticklabels()[::4]:
         label.set_visible(True)

   elif field[0:3] == 'ref':
      cbar.set_label('dBZ')

      for label in cbar.ax.xaxis.get_ticklabels():
         label.set_visible(True)

   elif field[2:6] == 'cape':
      cbar.set_label('J $\mathregular{kg^{-1}}$')

      for label in cbar.ax.xaxis.get_ticklabels():
         label.set_visible(True)

   elif field == 'max10mwind':
      cbar.set_label('kts')

      for label in cbar.ax.xaxis.get_ticklabels():
         label.set_visible(True)


   fig.suptitle(titlestr, size=18, weight='bold')

   fname = str.lower(model_str)+'_'+str.lower(field)+'_'+str.lower(domain)+'_'+cycle+'_f'+str(fhrs[0]-1).zfill(2)+'-f'+str(fhrs[-1]).zfill(2)

   plt.tight_layout()
   fig.subplots_adjust(top=0.96,bottom=0.04)
   plt.savefig(GRAPHX_DIR+'/'+fname+'.png',bbox_inches='tight')
   plt.close()



