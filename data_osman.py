import numpy as np
import csv

# Figure 1
cCl = 0.67976
Gamma_Li_0 = 459.74
CA0 = 0.64
with open('data/osman fig1 data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    A0_list_data = []
    Gamma_A_list_data = []
    cLi_list_data = []
    for row in datafile:
        A0_list_data += [float(row[0])]
        Gamma_A_list_data += [float(row[1])]
        cLi_list_data += [float(row[2])]
Gamma_Li_list_data = Gamma_Li_0*np.ones(len(cLi_list_data)) - cLi_list_data
cA_list_Data = CA0*(1-np.divide(Gamma_A_list_data, A0_list_data))
frac_A_list_data_fig1 = cA_list_Data/cCl
frac_ads_list_data_fig1 = Gamma_A_list_data/(Gamma_A_list_data+Gamma_Li_list_data)

# Figure 2
cCl = 0.6435
Gamma_Li_0 = 468.8737
CA0 = 0.62
with open('data/osman fig2 data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    A0_list_data = []
    Gamma_A_list_data = []
    cLi_list_data = []
    for row in datafile:
        A0_list_data += [float(row[0])]
        Gamma_A_list_data += [float(row[1])]
        cLi_list_data += [float(row[2])]
Gamma_Li_list_data = Gamma_Li_0*np.ones(len(cLi_list_data)) - cLi_list_data
cA_list_Data = CA0*(1-np.divide(Gamma_A_list_data, A0_list_data))
frac_A_list_data_fig2 = cA_list_Data/cCl
frac_ads_list_data_fig2 = Gamma_A_list_data/(Gamma_A_list_data+Gamma_Li_list_data)

# Figure 3
cCl = 0.59152
Gamma_Li_0 = 401.5072
CA0 = 0.56
with open('data/osman fig3 data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    A0_list_data = []
    Gamma_A_list_data = []
    cLi_list_data = []
    for row in datafile:
        A0_list_data += [float(row[0])]
        Gamma_A_list_data += [float(row[1])]
        cLi_list_data += [float(row[2])]
Gamma_Li_list_data = Gamma_Li_0*np.ones(len(cLi_list_data)) - cLi_list_data
cA_list_Data = CA0*(1-np.divide(Gamma_A_list_data, A0_list_data))
frac_A_list_data_fig3 = cA_list_Data/cCl
frac_ads_list_data_fig3 = Gamma_A_list_data/(Gamma_A_list_data+Gamma_Li_list_data)

# Figure 4
cCl = 0.66163
Gamma_Li_0 = 608.0727
CA0 = 0.67
with open('data/osman fig4 data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    A0_list_data = []
    Gamma_A_list_data = []
    cLi_list_data = []
    for row in datafile:
        A0_list_data += [float(row[0])]
        Gamma_A_list_data += [float(row[1])]
        cLi_list_data += [float(row[2])]
Gamma_Li_list_data = Gamma_Li_0*np.ones(len(cLi_list_data)) - cLi_list_data
cA_list_Data = CA0*(1-np.divide(Gamma_A_list_data, A0_list_data))
frac_A_list_data_fig4 = cA_list_Data/cCl
frac_ads_list_data_fig4 = Gamma_A_list_data/(Gamma_A_list_data+Gamma_Li_list_data)
