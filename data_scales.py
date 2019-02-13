import numpy as np
import csv

# Salt list
salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']

# Initiate dictionaries
c_data_fig1 = {}
zeta_data_fig1 = {}

# Extract data from CSVs
for salt in salts:
    with open('data/scales_fig1_data_'+salt+'.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        c = []
        zeta = []
        for row in datafile:
            c += [float(row[0])]
            zeta += [float(row[1])]
    c = np.array(c)
    zeta = np.array(zeta)
    c_data_fig1[salt] = c
    zeta_data_fig1[salt] = zeta
