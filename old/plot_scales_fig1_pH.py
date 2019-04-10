import matplotlib.pyplot as plt
import csv
import itertools
from scipy.optimize import least_squares
from solvers import *
from old.constants import *

# Parameters
salt_list = ['LiCl', 'NaCl', 'KCl', 'CsCl']
pH = 5.8
pKa = 5.6
C_1 = 10
C_2 = 2

marker = itertools.cycle(('o', 'v', 's', '8', 'p', '^', '<', '>'))
line = itertools.cycle(("-", "--", ":", "-."))

# Initialize lines and points
lines = []
points = []

# Extract data from csv and plot
data_dict = {}
for salt in salt_list:
    with open('data/scales_fig1_data_'+salt+'.csv') as csvfile:
        datafile = csv.reader(csvfile, delimiter=',', quotechar='|')
        conc_data = []
        zeta_data = []
        for row in datafile:
            conc_data += [float(row[0])]
            zeta_data += [float(row[1])]
        conc_data = np.array(conc_data)
        zeta_data = np.array(zeta_data)
    data_dict[salt] = np.vstack((conc_data, zeta_data))
    plt.scatter(np.log10(conc_data), zeta_data, marker=next(marker), label=salt+' (data)')

# Nonlinear regression
K_ads_dict = {}
X0 = np.array([3.1064756499971304])
regress = False
if regress:
    for salt in salt_list:
        def fun(X, c_data, zeta_data):
            K_ads = 10**X[0]
            res = np.zeros(len(c_data))
            i = 0
            for c in c_data:
                sol = solver_1plate_1cation_ads_pH(c, pH, K_ads, pKa, C_1=C_1, C_2=C_2)
                res[i] = sol.psi[2]*1000 - zeta_data[i]
                i += 1
            return res
        conc_data = data_dict[salt][0]
        zeta_data = data_dict[salt][1]
        res_lsq = least_squares(
            fun, X0,
            args=(conc_data, zeta_data),
            bounds=([-np.inf, np.inf]))
        K_ads_dict[salt] = res_lsq.x[0]
    print(res_lsq.x[0])

# Plot fit from model
K_ads_dict['LiCl'] = 10**-6.6510900712690795
K_ads_dict['NaCl'] = 10**-6.212985136745596
K_ads_dict['KCl'] = 10**2.5325569158774695
K_ads_dict['CsCl'] = 10**3.0689113645718678

for salt in salt_list:
    K_ads = K_ads_dict[salt]
    c_list = 10**np.linspace(-5, -1.5, 10)
    zeta_list = np.zeros(len(c_list))
    i = 0
    for c in c_list:
        sol = solver_1plate_1cation_ads_pH(c, pH, K_ads, pKa, C_1=C_1, C_2=C_2)
        zeta_list[i] = sol.psi[2]
        i += 1
    plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(line), label=salt+' (model)')
    font_size = 15
    plt.rc('font', size=font_size)
    plt.rc('axes', labelsize=font_size)
    plt.rc('xtick', labelsize=font_size)
    plt.rc('ytick', labelsize=font_size)
    plt.xlabel('log$_{10}(c)$ [mol/L]', fontsize=font_size)
    plt.ylabel('Zeta potential [mV]', fontsize=font_size)

plt.legend()
plt.show()
