import matplotlib.pyplot as plt
from scipy.constants import e, k, epsilon_0, N_A
from scipy.stats import linregress
from data_pashley import *


# Compute kappa
def compute_kappa(c, T=298, epsilon=79):
    rho = 1000*N_A*c
    return np.sqrt(2*rho*e**2/(epsilon*epsilon_0*k*T))


# Plot control
plot_fig1 = True
plot_fig2 = True
plot_fig3 = True
plot_fig4 = True

# Figure 1
if plot_fig1:
    plt.figure(1)
    pashley_fig1_c = []
    pashley_fig1_slope = []
    pashley_fig1_kappa = []
    start_index = [5, 6, 1, 5]
    i = 0
    for data in pashley_fig1_data:
        start = start_index[i]
        pashley_fig1_c += [data.c]
        pashley_fig1_kappa += [compute_kappa(data.c)]
        slope, intercept, r_value, p_value, std_err = linregress(data.D[start:], np.log(data.FR[start:]))
        pashley_fig1_slope += [slope]
        i += 1
        plt.scatter(data.D, data.FR, label='%.0e' % data.c)
        plt.semilogy(data.D, np.exp(intercept)*np.exp(slope*data.D))
    pashley_fig1_c = np.array(pashley_fig1_c)
    pashley_fig1_slope = np.array(pashley_fig1_slope)
    pashley_fig1_kappa = np.array(pashley_fig1_kappa)
    plt.title('Pashley fig 1')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N/m]')
    plt.legend(title='Concentration [M]')

# Figure 2
if plot_fig2:
    plt.figure(2)
    pashley_fig2_c = []
    pashley_fig2_slope = []
    pashley_fig2_kappa = []
    start_index = [0, 0, 0, 0]
    i = 0
    for data in pashley_fig2_data:
        start = start_index[i]
        pashley_fig2_c += [data.c]
        pashley_fig2_kappa += [compute_kappa(data.c)]
        slope, intercept, r_value, p_value, std_err = linregress(data.D[start:], np.log(data.FR[start:]))
        pashley_fig2_slope += [slope]
        i += 1
        plt.scatter(data.D, data.FR, label='%.0e' % data.c)
        plt.semilogy(data.D, np.exp(intercept)*np.exp(slope*data.D))
    pashley_fig2_c = np.array(pashley_fig2_c)
    pashley_fig2_slope = np.array(pashley_fig2_slope)
    pashley_fig2_kappa = np.array(pashley_fig2_kappa)
    plt.title('Pashley fig 2')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N/m]')
    plt.legend(title='Concentration [M]')

# Figure 3
if plot_fig3:
    plt.figure(3)
    pashley_fig3_c = []
    pashley_fig3_slope = []
    pashley_fig3_kappa = []
    start_index = [4, 8, 8]
    i = 0
    for data in pashley_fig3_data:
        start = start_index[i]
        pashley_fig3_c += [data.c]
        pashley_fig3_kappa += [compute_kappa(data.c)]
        slope, intercept, r_value, p_value, std_err = linregress(data.D[start:], np.log(data.FR[start:]))
        pashley_fig3_slope += [slope]
        i += 1
        plt.scatter(data.D, data.FR, label='%.0e' % data.c)
        plt.semilogy(data.D, np.exp(intercept)*np.exp(slope*data.D))
    pashley_fig3_c = np.array(pashley_fig3_c)
    pashley_fig3_slope = np.array(pashley_fig3_slope)
    pashley_fig3_kappa = np.array(pashley_fig3_kappa)
    plt.title('Pashley fig 3')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N/m]')
    plt.legend(title='Concentration [M]')

# Figure 4
if plot_fig4:
    plt.figure(4)
    pashley_fig4_c = []
    pashley_fig4_slope = []
    pashley_fig4_kappa = []
    start_index = [5, 6, 5]
    i = 0
    for data in pashley_fig4_data:
        start = start_index[i]
        pashley_fig4_c += [data.c]
        pashley_fig4_kappa += [compute_kappa(data.c)]
        slope, intercept, r_value, p_value, std_err = linregress(data.D[start:], np.log(data.FR[start:]))
        pashley_fig4_slope += [slope]
        i += 1
        plt.scatter(data.D, data.FR, label='%.0e' % data.c)
        plt.semilogy(data.D, np.exp(intercept)*np.exp(slope*data.D))
    pashley_fig4_c = np.array(pashley_fig4_c)
    pashley_fig4_slope = np.array(pashley_fig4_slope)
    pashley_fig4_kappa = np.array(pashley_fig4_kappa)
    plt.title('Pashley fig 4')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [$\mu$N/m]')
    plt.legend(title='Concentration [M]')

plt.show()
