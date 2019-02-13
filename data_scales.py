import numpy as np
import csv

# Figure 1
c_data_fig1 = {}
zeta_data_fig1 = {}
salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
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

# Figure 2
pH_data_fig2 = {}
zeta_data_fig2 = {}
with open('data/scales_fig2_data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    pH = []
    zeta = []
    for row in datafile:
        pH += [float(row[0])]
        zeta += [float(row[0])]
    pH_data_fig2 = np.array(pH)
    zeta_data_fig2 = np.array(zeta)
