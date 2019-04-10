import matplotlib.pyplot as plt
from old.data_scales import *
from old.class_Scales import *

# Set to true for plot
fig2 = True
fig10 = False
fig11 = False

# Figure 2
if fig2:
    plt.figure('plot_scales_fig2')
    plt.scatter(pH_data_fig2, zeta_data_fig2)
    fontsize = 14
    pH_list = np.linspace(2, 12)
    zeta_list = np.ones(len(pH_list))
    i = 0
    for pH in pH_list:
        sol = Solution_Scales(1e-3, 3.2, pKa=5.6, pH=pH)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(pH_list, zeta_list*1000)
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.title('Scales, Figure 2')
    plt.xlabel('pH', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.text(8, -20, 'p$K_a$ = %.1f\np*$K_M$ = %.1f' % (sol.pKa, sol.pKM))


# Figure 10
if fig10:
    from old.line_marker import *
    plt.figure('plot_scales_fig10')
    salts = ['CsCl', 'KCl']
    pKM = dict(CsCl=2.8, KCl=3.2)
    for salt in salts:
        c_data = c_data_fig1[salt]
        zeta_data = zeta_data_fig1[salt]
        plt.scatter(np.log10(c_data), zeta_data, marker=next(markers), label=salt)
        c_list = np.logspace(-4.5, -2)
        zeta_list = np.ones(len(c_list))
        i = 0
        for c in c_list:
            sol = Solution_Scales(c, pKM[salt])
            sol.solver_sigma()
            zeta_list[i] = sol.psi_d
            i += 1
        plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(lines))
    fontsize = 14
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.title('Scales, Figure 10')
    plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.text(-4.5, -55, 'pH=5.8')
    plt.legend()

# Figure 11
if fig11:
    from old.line_marker import *
    plt.figure('plot_scales_fig11')
    salts = ['NaCl', 'LiCl']
    for salt in salts:
        c_data = c_data_fig1[salt]
        zeta_data = zeta_data_fig1[salt]
        plt.scatter(np.log10(c_data), zeta_data, marker=next(markers), label=salt)
    c_list = np.logspace(-4, -2, 10)
    zeta_list = np.ones(len(c_list))
    i = 0
    for c in c_list:
        sol = Solution_Scales(c, 5.5)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    fontsize = 14
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.plot(np.log10(c_list), zeta_list*1000, color='k', linestyle=':')
    plt.title('Scales, Figure 11')
    plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.text(-4, -78, 'pH=5.8')
    plt.legend()

# Show all figures
plt.show()
