# transform geos-chem AOD output from utc time to local time

# libs
import scipy.io as sio
import numpy as np
import os

def Timeshift (path, path_, year, mon, nX, nY, nL, left_lon, xRes, day_ct):
    # dates
    if mon == '01' or mon == '03' or mon =='05' or mon == '07' or mon == '08' or mon == '10' or mon == '12':
        day = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

    if mon == '04' or mon == '06' or mon =='09' or mon == '11':
        day = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']

    if mon == '02':
        if int (year) % 4 == 0:
            day = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']
        else:
            day = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28']

    hour = ['0030', '0130', '0230', '0330', '0430', '0530', '0630', '0730', '0830', '0930', '1030', '1130', '1230', '1330', '1430', '1530', '1630', '1730', '1830', '1930', '2030', '2130', '2230', '2330']

    # composition
    comp = ['Sulfate', 'Nitrate', 'Ammonium', 'OA', 'BC', 'Dust', 'Seasalt', 'Total']
    tag = ['SpeciesConc', 'SpeciesConc', 'SpeciesConc', 'SpeciesConc', 'SpeciesConc', 'SpeciesConc', 'SpeciesConc', 'SpeciesConc']

    # time shift
    for iii in range (0, len (comp)):
        for j in range (0, day_ct):
            print (j)
            if j == 0:
                sourcefile = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile):
                    continue
                sourcedata = sio.loadmat (sourcefile)
                sourcedata = sourcedata [comp [iii]]

                sourcefile_plus = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j + 1] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile_plus):
                    continue
                sourcedata_plus = sio.loadmat (sourcefile_plus)
                sourcedata_plus = sourcedata_plus [comp [iii]]

                for k in range (12, 24):
                    outfile = path_ + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_' + hour [k] + 'l' + '.' + comp [iii] + '.AOD.mat'
                    data = np.zeros ((nX, nY, nL))

                    for ii in range (0, nX):
                        for jj in range (0, nY):
                            for ll in range (0, nL):
                                lon = left_lon + ii * xRes
                                utc = k + 0.5 - lon / 15.0
                        
                                if (utc >= 24.0):
                                    data [ii, jj, ll] = sourcedata_plus [int (utc - 24)] [ll] [jj] [ii]

                                else:
                                    data [ii, jj, ll] = sourcedata [int (utc)] [ll] [jj] [ii]
                
                    sio.savemat (outfile, {comp [iii]:data})

            elif j == day_ct - 1:
                sourcefile = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile):
                    continue
                sourcedata = sio.loadmat (sourcefile)
                sourcedata = sourcedata [comp [iii]]

                sourcefile_minus = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j - 1] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile_minus):
                    continue
                sourcedata_minus = sio.loadmat (sourcefile_minus)
                sourcedata_minus = sourcedata_minus [comp [iii]]

                for k in range (0, 12): 
                    outfile = path_ + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_' + hour [k] + 'l' + '.' + comp [iii] + '.AOD.mat'
                    data = np.zeros ((nX, nY, nL))
                
                    for ii in range (0, nX):
                        for jj in range (0, nY):
                            for ll in range (0, nL):
                                lon = left_lon + ii * xRes
                                utc = k + 0.5 - lon / 15.0

                                if (utc < 0.0):
                                    data [ii, jj, ll] = sourcedata_minus [int (utc + 24)] [ll] [jj] [ii]

                                else:
                                    data [ii, jj, ll] = sourcedata [int (utc)] [ll] [jj] [ii]
  
                    sio.savemat (outfile, {comp [iii]:data})                 

            else:
                #if (year + mon + day [j]) == '20170615':
                #    continue

                #if (year + mon + day [j - 1]) == '20170615':
                #    continue

                #if (year + mon + day [j + 1]) == '20170615':
                #    continue

                sourcefile = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile):
                    continue
                sourcedata = sio.loadmat (sourcefile)
                sourcedata = sourcedata [comp [iii]]

                sourcefile_plus = path + 'GEOSChem.'+ tag [iii] + '.' + year + mon + day [j + 1] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile_plus):
                    continue
                sourcedata_plus = sio.loadmat (sourcefile_plus)
                sourcedata_plus = sourcedata_plus [comp [iii]]

                sourcefile_minus = path + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j - 1] + '_0000z.nc4.AOD.' + comp [iii] + '.mat'
                if not os.path.exists (sourcefile_minus):
                    continue
                sourcedata_minus = sio.loadmat (sourcefile_minus)
                sourcedata_minus = sourcedata_minus [comp [iii]]

                for k in range (0, 24):
                    outfile = path_ + 'GEOSChem.' + tag [iii] + '.' + year + mon + day [j] + '_' + hour [k] + 'l' + '.' + comp [iii] + '.AOD.mat'
                    data = np.zeros ((nX, nY, nL))

                    for ii in range (0, nX):
                        for jj in range (0, nY):
                            for ll in range (0, nL):
                                lon = left_lon + ii * xRes
                                utc = k + 0.5 - lon / 15.0

                                if (utc >= 24.0): 
                                    data [ii, jj, ll] = sourcedata_plus [int (utc - 24)] [ll] [jj] [ii]
                
                                elif (utc < 0.0):
                                    data [ii, jj, ll] = sourcedata_minus [int (utc + 24)] [ll] [jj] [ii]
                                else:
                                    data [ii, jj, ll] = sourcedata [int (utc)] [ll] [jj] [ii]
                   
                    sio.savemat (outfile, {comp [iii]:data})

year = '2016'
mons = ['01']
day_ct = [31]
pathgc = '/storage1/fs1/rvmartin/Active/yanshun.li/GCC12.6.0/RunDir/na/geosfp/1.Alt/'
nX = 225
nY = 202
nL = 47
left_lon = -130.15625
xRes = 0.3125

for mm in range (0, len (mons)):
    pathin = pathgc + year + '/' + mons [mm] + '/AOD_' + str (day_ct [mm]) + 'days/'
    pathout = pathgc + year + '/' + mons [mm] + '/AOD_LT_' + str (day_ct [mm]) + 'days/'
    if not os.path.exists (pathout):
        os.mkdir (pathout)

    Timeshift (pathin, pathout, year, mons [mm], nX, nY, nL, left_lon, xRes, day_ct [mm])
