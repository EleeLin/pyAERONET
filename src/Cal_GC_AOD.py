# Calculate AOD from GEOS-Chem
import numpy as np
import os
from netCDF4 import Dataset
import scipy.io as sio
import glob
import math
from scipy.interpolate import interp1d

# Extract PM2.5 information from GEOS-Chem
def cal_gc_pm25 (pathin, pathmet, pathout, year, mon, nX, nY, day_ct):
    # Load optics
    OptTableDir = '/storage1/fs1/rvmartin/Active/GEOS-Chem-shared/ExtData/CHEM_INPUTS/FAST_JX/v2020-02/'

    AERName = ['so4', 'soot', 'org', 'ssa', 'ssc'] #[SO4/NH4/NIT, BC, OM/SOA, SALA, SALC]
    Dens = [1700, 1000, 1800, 2200, 2200] #unit in kg/m3, [SO4/NH4/NIT, BC, OM/SOA, SALA, SALC]
    RHbin = [0, 50, 70, 80, 90, 95, 99]
    #RHbin_SIA = [0, 35, 50, 70, 80, 90, 95, 99]
    Spec = [
    'SO4', 'NH4', 'NIT',\
    'BCPO', 'BCPI',\
    'OCPO', 'OCPI',\
    'SOAS', \
    'SALA',\
    'SALC'
    ]
    SpecGroup = np.concatenate([0 * np.ones(3), 1 * np.ones(2), 2 * np.ones(3), 3 * np.ones(1), 4 * np.ones(1)])
    SpecGroup = SpecGroup.astype(np.int16)
    #Hydrophilic = np.concatenate([True * np.ones(3), False * np.ones(1), True * np.ones(1), False * np.ones(1), True * np.ones(4)])
    #Hydrophilic = Hydrophilic.astype(np.int16)

    #Read scattering properties of Q, Reff at wavelength 550nm
    Q = np.zeros((len(AERName), len(RHbin)))
    Reff = np.zeros((len(AERName), len(RHbin)))

    for i in range(len(AERName)):
        fid = open(OptTableDir + AERName[i] + '.dat')
        lines = fid.read().splitlines()
        j = 0
        for line in lines:
            if line[0:5] == '  550':
                val = line.split()
                Q[i, j] = float(val[7])
                Reff[i, j] = float(val[9])
                j = j + 1
        fid.close()

    # For SIA: 0 growth at RH 35%
    #Q_SIA = np.concatenate([[Q[0, 0], Q[0, 0]], Q[0, 1:7]])
    #Reff_SIA = np.concatenate([[Reff[0, 0], Reff[0, 0]], Reff[0, 1:7]])

    #DST1 will be seperated into 4 different type with Mw = 29 * DustFrac (mass fraction)
    Dustspec = ['DST1', 'DST1', 'DST1', 'DST1', 'DST2', 'DST3', 'DST4']
    #DustFrac = [0.06, 0.12, 0.24, 0.58] # mass fraction
    DustFrac = [0.007, 0.0332, 0.2487, 0.7111] # mass fraction, updated accodring to L. Zhang 2013
    MwDust = 29
    DensDust = [2500, 2650] #Density for [DST1, DST2/DST3/DST4]

    # get QDust, ReffDust from dust.dat table
    QDust = np.zeros(len(Dustspec))
    ReffDust = np.zeros(len(Dustspec))

    fid = open(OptTableDir + 'dust.dat')
    lines = fid.read().splitlines()
    j = 0
    for line in lines:
        if line[0:5] == '  550':
            val = line.split()
            QDust[j] = float(val[7])
            ReffDust[j] = float(val[9])
            j = j + 1
    fid.close()

    # Search for GEOS-Chem diagnostic files
    if mon == '01' or mon == '03' or mon =='05' or mon == '07' or mon == '08' or mon == '10' or mon == '12':
        files = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
        ddd = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']

    if mon == '04' or mon == '06' or mon =='09' or mon == '11':
        files = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
        ddd = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']

    if mon == '02':
        if int (year) % 4 == 0:
            files = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ddd = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29']        
        else:
            files = ['', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']
            ddd = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28']

    for ff in range (0, len (ddd)):
        files [ff] = pathin + 'GEOSChem.SpeciesConc.' + year + mon + ddd [ff] + '_0000z.nc4'

    for dd in range (0, day_ct):
        filein = files [dd]
        print (filein)
        if (mon == '01' and ddd [dd] == '07'):
            continue

        filename_speciesconc = filein [len (pathin):len (filein)]
 
        fileout_Sulfate = pathout + filename_speciesconc + '.AOD.Sulfate.mat'  
        fileout_Nitrate = pathout + filename_speciesconc + '.AOD.Nitrate.mat'
        fileout_Ammonium = pathout + filename_speciesconc + '.AOD.Ammonium.mat'
        fileout_SIA = pathout + filename_speciesconc + '.AOD.SIA.mat'
        fileout_OA = pathout + filename_speciesconc + '.AOD.OA.mat'
        fileout_BC = pathout + filename_speciesconc + '.AOD.BC.mat' 
        fileout_Dust = pathout + filename_speciesconc + '.AOD.Dust.mat'
        fileout_Seasalt = pathout + filename_speciesconc + '.AOD.Seasalt.mat'
        fileout_Total = pathout + filename_speciesconc + '.AOD.Total.mat'
        fileout_OCPO = pathout + filename_speciesconc + '.AOD.OCPO.mat'
        fileout_BCPO = pathout + filename_speciesconc + '.AOD.BCPO.mat'

        filename_met = 'GEOSChem.StateMet.' + filename_speciesconc [21:35] + '.nc4'

        # retrieve data from GEOS-Chem
        gc_species_data = Dataset (filein, 'r')
        gc_met_data = Dataset (pathmet + filename_met, 'r')

        # SNA
        NH4 = np.array (gc_species_data.variables ['SpeciesConc_NH4'])
        NIT = np.array (gc_species_data.variables ['SpeciesConc_NIT'])
        SO4 = np.array (gc_species_data.variables ['SpeciesConc_SO4'])
        
        # BC
        BCPI = np.array (gc_species_data.variables ['SpeciesConc_BCPI'])
        BCPO = np.array (gc_species_data.variables ['SpeciesConc_BCPO'])

        # OC
        OCPI = np.array (gc_species_data.variables ['SpeciesConc_OCPI'])
        OCPO = np.array (gc_species_data.variables ['SpeciesConc_OCPO'])

        # DUST
        DST1 = np.array (gc_species_data.variables ['SpeciesConc_DST1'])
        DST2 = np.array (gc_species_data.variables ['SpeciesConc_DST2'])
        #DST3 = np.array (gc_species_data.variables ['SpeciesConc_DST3'])
        #DST4 = np.array (gc_species_data.variables ['SpeciesConc_DST4'])

        # SEASALT
        SALA = np.array (gc_species_data.variables ['SpeciesConc_SALA'])
        #SALC = np.array (gc_species_data.variables ['SpeciesConc_SALC'])

        # SOA
        SOAS = np.array (gc_species_data.variables ['SpeciesConc_SOAS'])

        # Met
        PMID = np.array (gc_met_data.variables ['Met_PMIDDRY'])
        T = np.array (gc_met_data.variables ['Met_T'])
        SPHU = np.array (gc_met_data.variables ['Met_SPHU'])
        BxH = np.array (gc_met_data.variables ['Met_BXHEIGHT'])
       
        RH = np.zeros ((24, 47, nY, nX))
        
        for i in range (0, 24):
            for j in range (0, 47):
                for k in range (0, nY):
                    for l in range (0, nX):
                        RH [i, j, k, l] = 0.263 * PMID [i, j, k, l] * 100.0 * SPHU [i, j, k, l] / 1000.0 / np.exp (17.625 * (T [i, j, k, l] - 273.16) / (T [i, j, k, l] - 29.65))
       
        #print (np.max (RH))
        #print (np.min (RH))
 
        RH [RH > 90] = 90

        # convert mixing ratio to kg/m3
        ppb_ugm3  = 1e6 / 8.314 * 100.0 * PMID / T * 1.0e-9

        # SNA
        NH4 = NH4 * ppb_ugm3 * 18.0
        NIT = NIT * ppb_ugm3 * 62.0
        SO4 = SO4 * ppb_ugm3 * 96.0

        # BC
        BCPI = BCPI * ppb_ugm3 * 12.01
        BCPO = BCPO* ppb_ugm3 * 12.01

        # OC
        OCPI = OCPI * ppb_ugm3 * 12.01
        OCPO = OCPO * ppb_ugm3 * 12.01
               
        # DUST
        DST1 = DST1 * ppb_ugm3 * 29.0
        DST2 = DST2 * ppb_ugm3 * 29.0
        #DST3 = DST3 * ppb_ugm3 * 29.0
        #DST4 = DST4 * ppb_ugm3 * 29.0

        # SEASALT
        SALA = SALA * ppb_ugm3 * 31.4
        #SALC = SALC * ppb_ugm3 * 31.4

        # SOA
        SOAS = SOAS * ppb_ugm3 * 150.0

        # Non dust AOD
        NH4_Rdry = Reff [0, 0] * 1e-6 #unit is m
        NH4_Qdry = Q [0, 0]
        NH4_fQ = interp1d (RHbin, Q [0, :], kind='linear', fill_value='extrapolate')
        NH4_fR = interp1d (RHbin, Reff [0, :], kind='linear', fill_value='extrapolate')
        NH4_tQ = NH4_fQ (RH)
        NH4_tReff = NH4_fR (RH) * 1.0e-6 # unit is m
        NH4_AOD = 0.75 * NH4_tQ * NH4 * BxH / (Dens [0] * NH4_Rdry) * ((NH4_tReff / NH4_Rdry) ** 2)

        SO4_Rdry = Reff [0, 0] * 1e-6 #unit is m
        SO4_Qdry = Q [0, 0]
        SO4_fQ = interp1d (RHbin, Q [0, :], kind='linear', fill_value='extrapolate')
        SO4_fR = interp1d (RHbin, Reff [0, :], kind='linear', fill_value='extrapolate')
        SO4_tQ = SO4_fQ (RH)
        SO4_tReff = SO4_fR (RH) * 1.0e-6 # unit is m
        SO4_AOD = 0.75 * SO4_tQ * SO4 * BxH / (Dens [0] * SO4_Rdry) * ((SO4_tReff / SO4_Rdry) ** 2)

        NIT_Rdry = Reff [0, 0] * 1e-6 #unit is m
        NIT_Qdry = Q [0, 0]
        NIT_fQ = interp1d (RHbin, Q [0, :], kind='linear', fill_value='extrapolate')
        NIT_fR = interp1d (RHbin, Reff [0, :], kind='linear', fill_value='extrapolate')
        NIT_tQ = NIT_fQ (RH)
        NIT_tReff = NIT_fR (RH) * 1.0e-6 # unit is m
        NIT_AOD = 0.75 * NIT_tQ * NIT * BxH / (Dens [0] * NIT_Rdry) * ((NIT_tReff / NIT_Rdry) ** 2)

        BCPI_Rdry = Reff [1, 0] * 1e-6 #unit is m
        BCPI_Qdry = Q [1, 0]
        BCPI_fQ = interp1d (RHbin, Q [1, :], kind='linear', fill_value='extrapolate')
        BCPI_fR = interp1d (RHbin, Reff [1, :], kind='linear', fill_value='extrapolate')
        BCPI_tQ = BCPI_fQ (RH)
        BCPI_tReff = BCPI_fR (RH) * 1.0e-6 # unit is m
        BCPI_AOD = 0.75 * BCPI_tQ * BCPI * BxH / (Dens [1] * BCPI_Rdry) * ((BCPI_tReff / BCPI_Rdry) ** 2)

        BCPO_Rdry = Reff [1, 0] * 1e-6 #unit is m
        BCPO_Qdry = Q [1, 0]
        BCPO_AOD = 0.75 * BCPO_Qdry * BCPO * BxH / (Dens [1] * BCPO_Rdry)

        OCPO_Rdry = Reff [2, 0] * 1e-6 #unit is m
        OCPO_Qdry = Q [2, 0]
        OCPO_AOD = 0.75 * OCPO_Qdry * OCPO * 2.1 * BxH / (Dens [2] * OCPO_Rdry)

        OCPI_Rdry = Reff [2, 0] * 1e-6 #unit is m
        OCPI_Qdry = Q [2, 0]
        OCPI_fQ = interp1d (RHbin, Q [2, :], kind='linear', fill_value='extrapolate')
        OCPI_fR = interp1d (RHbin, Reff [2, :], kind='linear', fill_value='extrapolate')
        OCPI_tQ = OCPI_fQ (RH)
        OCPI_tReff = OCPI_fR (RH) * 1.0e-6 # unit is m
        OCPI_AOD = 0.75 * OCPI_tQ * OCPI * 2.1 * BxH / (Dens [2] * OCPI_Rdry) * ((OCPI_tReff / OCPI_Rdry) ** 2)

        SOAS_Rdry = Reff [2, 0] * 1e-6 #unit is m
        SOAS_Qdry = Q [2, 0]
        SOAS_fQ = interp1d (RHbin, Q [2, :], kind='linear', fill_value='extrapolate')
        SOAS_fR = interp1d (RHbin, Reff [2, :], kind='linear', fill_value='extrapolate')
        SOAS_tQ = SOAS_fQ (RH)
        SOAS_tReff = SOAS_fR (RH) * 1.0e-6 # unit is m
        SOAS_AOD = 0.75 * SOAS_tQ * SOAS * BxH / (Dens [2] * SOAS_Rdry) * ((SOAS_tReff / SOAS_Rdry) ** 2)

        SALA_Rdry = Reff [3, 0] * 1e-6 #unit is m
        SALA_Qdry = Q [3, 0]
        SALA_fQ = interp1d (RHbin, Q [3, :], kind='linear', fill_value='extrapolate')
        SALA_fR = interp1d (RHbin, Reff [3, :], kind='linear', fill_value='extrapolate')
        SALA_tQ = SALA_fQ (RH)
        SALA_tReff = SALA_fR (RH) * 1.0e-6 # unit is m
        SALA_AOD = 0.75 * SALA_tQ * SALA * BxH / (Dens [3] * SALA_Rdry) * ((SALA_tReff / SALA_Rdry) ** 2)        

        #SALC_Rdry = Reff [4, 0] * 1e-6 #unit is m
        #SALC_Qdry = Q [4, 0]
        #SALC_fQ = interp1d (RHbin, Q [4, :], kind='linear', fill_value='extrapolate')
        #SALC_fR = interp1d (RHbin, Reff [4, :], kind='linear', fill_value='extrapolate')
        #SALC_tQ = SALC_fQ (RH)
        #SALC_tReff = SALC_fR (RH) * 1.0e-6 # unit is m
        #SALC_AOD = 0.75 * SALC_tQ * SALC * BxH / (Dens [4] * SALC_Rdry) * ((SALC_tReff / SALC_Rdry) ** 2)     

        # Dust AOD
        DST1_frac1_AOD = 0.75 * QDust [0] * DST1 * DustFrac [0] * BxH / (DensDust[0] * ReffDust [0] * 1.0e-6)
        DST1_frac2_AOD = 0.75 * QDust [1] * DST1 * DustFrac [1] * BxH / (DensDust[0] * ReffDust [1] * 1.0e-6)
        DST1_frac3_AOD = 0.75 * QDust [2] * DST1 * DustFrac [2] * BxH / (DensDust[0] * ReffDust [2] * 1.0e-6)
        DST1_frac4_AOD = 0.75 * QDust [3] * DST1 * DustFrac [3] * BxH / (DensDust[0] * ReffDust [3] * 1.0e-6)
        DST2_AOD = 0.75 * QDust [4] * DST2 * BxH / (DensDust[1] * ReffDust [4] * 1.0e-6)
        #DST3_AOD = 0.75 * QDust [5] * DST3 * BxH / (DensDust[1] * ReffDust [5] * 1.0e-6)
        #DST4_AOD = 0.75 * QDust [6] * DST4 * BxH / (DensDust[1] * ReffDust [6] * 1.0e-6)

        SIA_AOD = SO4_AOD + NIT_AOD + NH4_AOD
        BC_AOD = BCPI_AOD + BCPO_AOD
        OM_AOD = OCPI_AOD + OCPO_AOD + SOAS_AOD
        SS_AOD = SALA_AOD# + SALC_AOD
        DUST_AOD = DST1_frac1_AOD + DST1_frac2_AOD + DST1_frac3_AOD + DST1_frac4_AOD + DST2_AOD# + DST3_AOD + DST4_AOD
        TOTAL_AOD = SIA_AOD + BC_AOD + OM_AOD + SS_AOD + DUST_AOD

        sio.savemat (fileout_Sulfate, {'Sulfate':SO4_AOD})
        sio.savemat (fileout_Nitrate, {'Nitrate':NIT_AOD})
        sio.savemat (fileout_Ammonium, {'Ammonium':NH4_AOD})
        sio.savemat (fileout_OA, {'OA':OM_AOD})
        sio.savemat (fileout_BC, {'BC':BC_AOD})
        sio.savemat (fileout_Dust, {'Dust':DUST_AOD})
        sio.savemat (fileout_Seasalt, {'Seasalt':SS_AOD})
        sio.savemat (fileout_Total, {'Total':TOTAL_AOD})
        sio.savemat (fileout_OCPO, {'OCPO':OCPO_AOD})
        sio.savemat (fileout_BCPO, {'BCPO':BCPO_AOD})

# main
#mons = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
mons = ['01']
#day_ct = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
day_ct = [31]
pathin = '/storage1/fs1/rvmartin/Active/yanshun.li/GCC12.6.0/RunDir/na/geosfp/1.Alt/2016/'
pathmet = '/storage1/fs1/rvmartin/Active/yanshun.li/GCC12.6.0/RunDir/na/geosfp/9.Met/2016/'
nX = 225
nY = 202

for mm in range (0, len (mons)):
    pathout = pathin + mons [mm] + '/AOD_' + str (day_ct [mm]) + 'days/'
    if not os.path.exists (pathout):
        os.mkdir (pathout)
    cal_gc_pm25 (pathin + mons [mm] + '/geosfp_025x03125_tropchem_na/OutputDir/', pathmet + mons [mm] + '/geosfp_025x03125_tropchem_na/OutputDir/', pathout, '2016', mons [mm], nX, nY, day_ct [mm]) 
