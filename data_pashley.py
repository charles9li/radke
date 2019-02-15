import numpy as np
import csv


# Class for storing data
class PashleyData:
    def __init__(self, c, D, FR):
        self.c = c
        self.D = D
        self.FR = FR


# Figure 1
c_list = ['6e-2', '1e-2', '1e-3', '1e-4']
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
