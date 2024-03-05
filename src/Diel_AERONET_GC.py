# libs
import csv
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 25})
import numpy as np
from matplotlib.ticker import FormatStrFormatter

# regions
regions = ['West', 'Central', 'East']

# seasons
seasons = ['DJF', 'MAM', 'JJA', 'SON']

for rr in range (0, len (regions)):
    for ss in range (0, len (seasons)):
        fig, ax = plt.subplots (figsize = (10.0, 8.0))
        plt.subplots_adjust (left = 0.20)

        # diel norm data
        filename = './diel_' + regions [rr] + '_' + seasons [ss] + '.csv'
        reader = csv.reader (open (filename, "rt"), delimiter = ",")
        diel_data = list (reader)

        aeronet_diel = diel_data [0]
        goesr_diel = diel_data [1]
        ct_diel = diel_data [2]
        indicator = np.zeros (24)
        for i in range (0, len (aeronet_diel)):
            aeronet_diel [i] = float (aeronet_diel [i])
            goesr_diel [i] = float (goesr_diel [i])
            ct_diel [i] = float (ct_diel [i])
            if ct_diel [i] < 100.0:
                aeronet_diel [i] = np.nan
                goesr_diel [i] = np.nan
                indicator [i] = 1

        aeronet_diel_ = []
        goesr_diel_ = []
        for iii in range (0, 24):
            if np.isnan (aeronet_diel [iii]) == False and np.isnan (goesr_diel [iii]) == False:
                aeronet_diel_.append (aeronet_diel [iii])
                goesr_diel_.append (goesr_diel [iii])

        print (aeronet_diel_)
        print (goesr_diel_)
        RMSD = math.sqrt (np.square (np.subtract (aeronet_diel_, goesr_diel_)).mean())
        RMSD_str = f"{RMSD:.3f}"

        outfile = './diel_' + regions [rr] + '_' + seasons [ss] + '.jpg'

        #hour = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        hour = []
        x_data = []
        for iii in range (0, 24):
            if indicator [iii] == 0:
                x_data.append (iii)
                hour.append (iii)

        #plt.xticks (np.arange (min (hour) - 1, max (hour) + 1, 4.0))
        #plt.xticks ([0, 4, 8, 12, 16, 20, 23])
        plt.xticks (x_data)
        plt.plot (hour, aeronet_diel_, color = 'black', linewidth = 5.0)
        plt.plot (hour, goesr_diel_, color = 'red', linewidth = 5.0)
        plt.xlabel ('Local Solar Time (hr)')
        plt.ylabel ('AOD')
        plt.gca ().legend (('AERONET', 'GEOS-Chem (RMSD=' + RMSD_str + ')'), frameon = False)
        #plt.xlim ((-1, 24))
        plt.ylim ((0.0, 0.15))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.savefig (outfile)
        plt.clf ()
