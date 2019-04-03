import numpy as np
import csv


# Class for storing data
class ClaessonData:
    def __init__(self, wash, press):
        self.wash = wash
        self.press = press


# Figure 2
with open('data/claesson_fig2_data.csv') as csvfile:
    datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
    M_list = []
    S_list = []
    for row in datafile:
        M_list += [float(row[0])]
        S_list += [float(row[1])]
wash = np.vstack((np.array(M_list[0:9]), np.array(S_list[0:9])))
press = np.vstack((np.array(M_list[9:-1]), np.array(S_list[9:-1])))
claesson_fig2_data = ClaessonData(wash, press)
