import matplotlib.pyplot as plt
from scipy.constants import e, k, epsilon_0, N_A
from scipy.stats import linregress
from data_pashley import *


# Compute kappa
def compute_kappa(c, T=298, epsilon=79):
    rho = 1000*N_A*c
    return np.sqrt(2*rho*e**2/(epsilon*epsilon_0*k*T))


# Figure 1
plt.figure(1)
pashley_fig1_c = []
pashley_fig1_slope = []
pashley_fig1_kappa = []
start_index = [5, 0, 0, 5]
i = 0
for data in pashley_fig1_data:
    pashley_fig1_c += [data.c]
    pashley_fig1_kappa += [compute_kappa(data.c)]
    slope, intercept, r_value, p_value, std_err = linregress(data.D[i:], np.log(data.FR[i:]))
    pashley_fig1_slope += [slope]
    i += 1
    plt.scatter(data.D, data.FR)
    plt.semilogy(data.D, np.exp(intercept)*np.exp(slope*data.D))
pashley_fig1_c = np.array(pashley_fig1_c)
pashley_fig1_slope = np.array(pashley_fig1_slope)
pashley_fig1_kappa = np.array(pashley_fig1_kappa)

plt.show()
