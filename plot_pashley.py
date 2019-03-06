import matplotlib.pyplot as plt
from scipy.constants import e, k
from data_pashley import *
from class_Pashley import *

# Plot control
plot_fig1 = False
plot_fig2 = True
plot_fig3 = True

# Figure 1
if plot_fig1:
    # psi_d_list_fig1 = np.ones(len(pashley_fig1_data))
    # for data in pashley_fig1_data:
    #     c = data.c
    data = pashley_fig1_data[0]
    c = data.c
    D = data.D
    FR = data.FR
    plt.scatter(D, FR*1e-3)
    sol = PashleySolution(c, -.095, np.linspace(D[0], D[-1])*1e-9)
    sol.kappa = 1/21.08e-9
    sol.compute_FR()
    plt.semilogy(sol.D_list*1e9, sol.FR_list*1e3)
    print(c)

# Figure 2
if plot_fig2:
    plt.figure('2')
    c_list = []
    zeta_list = np.array([116.8, 120, 145.4, 70])/-1000
    kappa_list = 1e9/np.array([21.13, 16.24, 7.65, 2.66])
    for i in range(len(pashley_fig2_data)):
        data = pashley_fig2_data[i]
        c = data.c
        c_list += [float(c)]
        D = data.D
        FR = data.FR
        plt.scatter(D, FR*1e-3, label='%.1e' % c)
        zeta = zeta_list[i]
        sol = PashleySolution(c, zeta, np.linspace(D[0], D[-1])*1e-9)
        sol.kappa = kappa_list[i]
        sol.compute_FR()
        plt.semilogy(sol.D_list*1e9, sol.FR_list*1e3)
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.title('Pashley Fig 2 (NaCl)')
    plt.legend(title='concentration [M]')

    plt.figure('2_zeta')
    ax = plt.gca()
    ax.set_xscale('log')
    ax.scatter(c_dict['NaCl'], zeta_dict['NaCl'], label='Pashley, fig 2')
    ax.scatter(c_list, -zeta_list*1000, label='fits from data')
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel('bulk ion conc [mol/L]')
    ax.set_ylabel('surface potential [mV]')
    ax.set_title('NaCl')
    ax.legend()

# Figure 3
if plot_fig3:
    plt.figure('3')
    c_list = []
    zeta_list = np.array([85, 100, 90])/-1000
    kappa_list = 1e9/np.array([24.12, 14.7, 9.01])
    for i in range(len(pashley_fig3_data)):
        data = pashley_fig3_data[i]
        c = data.c
        c_list += [float(c)]
        D = data.D
        FR = data.FR
        plt.scatter(D, FR*1e-3, label='%.1e' % c)
        zeta = zeta_list[i]
        sol = PashleySolution(c, zeta, np.linspace(D[0], D[-1])*1e-9)
        sol.kappa = kappa_list[i]
        sol.compute_FR()
        plt.semilogy(sol.D_list*1e9, sol.FR_list*1e3)
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.title('Pashley Fig 3 (KCl)')
    plt.legend(title='concentration [M]')

    plt.figure('3_zeta')
    ax = plt.gca()
    ax.set_xscale('log')
    ax.scatter(c_dict['KCl'], zeta_dict['KCl'], label='Pashley, fig 3')
    ax.scatter(c_list, -zeta_list*1000, label='fits from data')
    ax.set_xlim(ax.get_xlim()[::-1])
    ax.set_xlabel('bulk ion conc [mol/L]')
    ax.set_ylabel('surface potential [mV]')
    ax.set_title('KCl')
    ax.legend()

plt.show()
