import numpy as np
import h5py
import csv
import math
import scipy.io as sio
import os.path

def readAERONET(year,month,hours,**kwargs):

    #hours in UTC
    file = '/storage1/fs1/rvmartin/Active/yanshun.li/Progess_reports/2022Aug19/GoogleEarthAERONET.20220208-AERONET-Hourly.mat'
    #hrAERdata=loadmat(file)

    hrAERdata = h5py.File(file)

    # station=np.array(hrAERdata['station'])
    # print(hrAERdata[station[0][0]])
    # exit()
    coord=np.array(hrAERdata['coord']) #station coordinates

    lats=coord[0,:]
    lons=coord[1,:]
    #lons[lons<0]=lons[lons<0]+360.

    #AWLNS=hrAERdata['AWLNS'] #all wavelengths provided directly by AERONET
    WLNS=np.array(hrAERdata['WLNS']) #satellite-relevant wavelengths calculated from observed AOD at AWVLNSdates: date definition dimension
    #adata=hrAERdata['adata'] #average hourly AERONET AOD, with dimensions [station,dates,awvlns,hour]
    #adatac=hrAERdata['adatac'] #number of observations included in adata
    data=np.array(hrAERdata['data']) #average hourly AERONET AOD calculated at WVLNS using angstrom exponent determined from adata, with dimensions [station,dates,wvlns,hour]
    #datac=hrAERdata['datac'] #number of observations included in data
    dates=np.array(hrAERdata['dates']).flatten()
    hrAERdata.close()

    years=np.floor(dates/10000).astype(int)
    months=np.floor((dates-years*10000)/100).astype(int)

    minyear=np.min(year)
    maxyear=np.max(year)

    minmonth=np.min(month)
    maxmonth=np.max(month)


    #confine the time period
    tinds=np.array(((years>=minyear)&(years<=maxyear)&(months>=minmonth)&(months<=maxmonth)).nonzero()).flatten()

    tdata=(data[:,0,tinds,:])[hours,:,:] #hours,wavelength,dates,sites
    tdates=dates[tinds]

    #confine the lat lon range
    if 'rglatlon' in kwargs:
        rglatlon=kwargs['rglatlon']
        minlon=rglatlon[0]
        maxlon=rglatlon[1]
        minlat=rglatlon[2]
        maxlat=rglatlon[3]
        siteinds=np.array(((lats>=minlat)&(lats<=maxlat)&(lons>=minlon)&(lons<=maxlon)).nonzero()).flatten()

        tdata=tdata[:,:,siteinds]
        # station=station[siteinds]
        lats=lats[siteinds]
        lons=lons[siteinds]

    return [lats,lons,tdates,tdata]

year = [2016]
month = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
hours = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

AERONET_data = readAERONET (year, month, hours)

hours_str = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23']
dayct_str = ['31', '29', '31', '30', '31', '30', '31', '31', '30', '31', '30', '31']
mon_str = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
gc_path = '/storage1/fs1/rvmartin/Active/yanshun.li/GCC12.6.0/RunDir/na/geosfp/1.Alt/'
matched_data = []

for ii in range (0, 24):
    for jj in range (0, 366):
        date = str (int (AERONET_data [2] [jj]))
        year = date [0:4]
        mon = date [4:6]
        day = date [6:8]
        hour = hours_str [ii]
        gc_file = gc_path + year + '/' + mon + '/AOD_LT_' + dayct_str [mon_str.index (mon)] + 'days/GEOSChem.SpeciesConc.' + date + '_' + hour + '30l.Total.AOD.mat'

        if (os.path.exists (gc_file) == False):
            continue

        gc_aod = sio.loadmat (gc_file)
        gc_aod = gc_aod ['Total']
        gc_aod_ = np.zeros ((225, 202))
        for mm in range (0, 225):
            for nn in range (0, 202):
                for pp in range (0, 47):
                    gc_aod_ [mm, nn] += gc_aod [mm, nn, pp] 
        
        print (date + '_' + hour)

        for kk in range (0, 1303):
            lat = float (AERONET_data [0] [kk])
            lon = float (AERONET_data [1] [kk])
            
            if (math.isnan (lon) == True or math.isnan (lat) == True):
                continue

            nX = int ((lon + 130.15625)  / 0.3125)
            nY = int ((lat - 9.625) / 0.25)

            if (nX > 224 or nX < 0 or nY > 201 or nY < 0):
                continue

            aeronet_aod_value = float (AERONET_data [3] [ii, jj, kk])
            gc_aod_value = gc_aod_ [nX, nY]

            if (math.isnan (aeronet_aod_value) == False and math.isnan (gc_aod_value) == False):
                if (aeronet_aod_value > 0.0 and gc_aod_value > 0.0):
                    matched_data.append ([date, hour, lon, lat, aeronet_aod_value, gc_aod_value])

outfile = open ('./MatchedAOD.csv', 'w')
with outfile:
    writer = csv.writer (outfile)
    for i in range (0, len (matched_data)):
        writer.writerow (matched_data [i])
