import numpy as np
import csv


# Class for storing data
class PashleyData:
    def __init__(self, c, D, FR):
        self.c = c
        self.D = D
        self.FR = FR


# Figure 1
c_list = ['1e-4', '1e-3', '1e-2', '6e-2']
pashley_fig1_data = []
for c in c_list:
    with open('data/pashley_fig1_data_' + c + '.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        D = []
        FR = []
        for row in datafile:
            D += [float(row[0])]
            FR += [float(row[1])]
    D = np.array(D)
    FR = np.array(FR)
    pashley_fig1_data += [PashleyData(float(c), D, FR)]

# Figure 2
c_list = ['4e-5', '1e-4', '1e-3', '1e-2']
pashley_fig2_data = []
for c in c_list:
    with open('data/pashley_fig2_data_' + c + '.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        D = []
        FR = []
        for row in datafile:
            D += [float(row[0])]
            FR += [float(row[1])]
    D = np.array(D)
    FR = np.array(FR)
    pashley_fig2_data += [PashleyData(float(c), D, FR)]

# Figure 3
c_list = ['4e-5', '3e-4', '1e-3']
pashley_fig3_data = []
for c in c_list:
    with open('data/pashley_fig3_data_' + c + '.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        D = []
        FR = []
        for row in datafile:
            D += [float(row[0])]
            FR += [float(row[1])]
    D = np.array(D)
    FR = np.array(FR)
    pashley_fig3_data += [PashleyData(float(c), D, FR)]

# Figure 4
c_list = ['4e-5', '1e-4', '1e-3']
pashley_fig4_data = []
for c in c_list:
    with open('data/pashley_fig4_data_' + c + '.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        D = []
        FR = []
        for row in datafile:
            D += [float(row[0])]
            FR += [float(row[1])]
    D = np.array(D)
    FR = np.array(FR)
    pashley_fig4_data += [PashleyData(float(c), D, FR)]

# Figure 5
c_dict = {}
zeta_dict = {}
salt_list = ['LiCl', 'NaCl', 'KCl', 'CsCl']
for salt in salt_list:
    with open('data/pashley_fig5_data_'+salt+'.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        c = []
        zeta = []
        for row in datafile:
            c += [float(row[0])]
            zeta += [float(row[1])]
    c_dict[salt] = np.array(c)
    zeta_dict[salt] = np.array(zeta)
