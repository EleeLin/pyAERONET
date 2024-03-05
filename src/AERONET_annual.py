import numpy as np
import h5py
import csv
import math

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
Annual_AOD = np.zeros (len (AERONET_data [0]))
Annual_ct = np.zeros (len (AERONET_data [1]))
for ii in range (0, 24):
    for jj in range (0, 366):
        for kk in range (0, 1303):
            aod_value = float (AERONET_data [3] [ii, jj, kk])
            if (math.isnan (aod_value) == False): 
                Annual_AOD [kk] += aod_value
                Annual_ct [kk] += 1.0

for i in range (0, len (AERONET_data [0])):
    if Annual_ct [i] != 0.0:
        Annual_AOD [i] /= Annual_ct [i]

outfile = open ('./AERONET_annual.csv', 'w')
with outfile:
    writer = csv.writer (outfile)
    for i in range (0, len (AERONET_data [0])):
        if Annual_AOD [i] != 0.0:
            writer.writerow ([AERONET_data [1] [i], AERONET_data [0] [i], Annual_AOD [i], Annual_ct [i] / 366.0 / 24.0])
