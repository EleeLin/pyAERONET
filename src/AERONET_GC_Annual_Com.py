import csv
import numpy as np

# load hourly comparison file
hourly_com_file = './MatchedAOD.csv'
reader = csv.reader (open (hourly_com_file, "rt"), delimiter = ",")
hourly_com_data = list (reader)

# get sites
site_pos = []
for i in range (0, len (hourly_com_data)):
    lon = float (hourly_com_data [i] [2])
    lat = float (hourly_com_data [i] [3])

    if [lon, lat] not in site_pos:
        site_pos.append ([lon, lat])

aeronet_aod = np.zeros (len (site_pos))
gc_aod = np.zeros (len (site_pos))
ct_aod = np.zeros (len (site_pos))

for i in range (0, len (hourly_com_data)):
    lon = float (hourly_com_data [i] [2])
    lat = float (hourly_com_data [i] [3])

    index = site_pos.index ([lon, lat])
    aeronet_aod [index] += float (hourly_com_data [i] [4])
    gc_aod [index] += float (hourly_com_data [i] [5])
    ct_aod [index] += 1.0

for i in range (0, len (site_pos)):
    if ct_aod [i] > 0.0:
        aeronet_aod [i] /= ct_aod [i]
        gc_aod [i] /= ct_aod [i]

outfile = open ('./AERONET_GC_Annual_Com_2016.csv', 'w')
with outfile:
    writer = csv.writer (outfile)
    for i in range (0, len (site_pos)):
        if float (site_pos [i] [0]) < -125.0 or float (site_pos [i] [0]) > -66.0 or float (site_pos [i] [1]) > 49.0 or float (site_pos [i] [1]) < 25.0:
            continue
        else: 
            writer.writerow ([site_pos [i] [0], site_pos [i] [1], aeronet_aod [i], gc_aod [i], ct_aod [i] / 366.0 / 24.0])
