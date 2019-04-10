import numpy as np
import csv
import matplotlib.pyplot as plt


# Class for storing data
class IsraelachviliData:
    def __init__(self, c, D, FR):
        self.c = c
        self.D = D
        self.FR = FR


# Figure 3
c_list = ['1e-4', '1e-3', '1e-2', '1e-1']
israelachvili_fig3_data = []
for c in c_list:
    with open('data/israelachvili_fig3_data_' + c + '.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        D = []
        FR = []
        for row in datafile:
            D += [float(row[0])]
            FR += [float(row[1])]
    D = np.array(D)
    FR = np.array(FR)
    israelachvili_fig3_data += [IsraelachviliData(float(c), D, FR)]
