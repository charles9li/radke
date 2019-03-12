import matplotlib.pyplot as plt
from data_scales import *
from data_osman import *
from line_marker import *
from radii import *
from solver import *


##############
# PARAMETERS #
##############

# K_ads
K_Li = 10**0.5
K_Na = 10**2
K_K = 10**2.8
K_Cs = 10**3
K_ads_list = [K_Li, K_Na, K_K, K_Cs]

# pH
pH = 5.8

# pKa
pKa = 5.3

# Permittivity
eps1 = 6*epsilon_0
eps2 = 30*epsilon_0
eps_bulk = 79*epsilon_0

# C1
C1_Li = eps1/R_Li
C1_Na = eps1/R_Na
C1_K = eps1/R_K
C1_Cs = eps1/R_Cs
C1_list = [C1_Li, C1_Na, C1_K, C1_Cs]

# C2
C2_Li = eps2/(2*R_Li_hyd)
C2_Na = eps2/(2*R_Na_hyd)
C2_K = eps2/(2*R_K_hyd)
C2_Cs = eps2/(2*R_Cs_hyd)
C2_list = [C2_Li, C2_Na, C2_K, C2_Cs]


#########
# PLOTS #
#########

# Scales fig 1
plot_scales_fig1 = True
if plot_scales_fig1:
    # plt.figure('Scales fig 1')
    i = 0
    # salts = ['LiCl', 'NaCl']
    for salt in salts:
        plt.figure('Scales fig 1 ' + salt)
        K_ads = K_ads_list[i]
        C1 = C1_list[i]
        C2 = C2_list[i]
        color = next(colors)
        c_list = 10**np.linspace(-5, -1, 20)
        zeta_array = np.zeros((3, len(c_list)))
        j = 0
        for c in c_list:
            k = 0
            for pH in np.array([5.5, 5.8, 6.1]):
                c_list1 = np.array([c+10**-pH, c, 10**-pH])
                K_list = np.array([K_ads, 10**pKa])
                z_list = np.array([-1, 1, 1])
                v_list = np.array([False, True, True])
                sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                                      pH=pH, pH_effect=False, C_1=C1, C_2=C2)
                sol.solver_sigma()
                zeta_array[k, j] = sol.psi_d
                k += 1
            j += 1
        plt.plot(np.log10(c_list), zeta_array[1]*1000, color='k', label='model')
        plt.fill_between(np.log10(c_list), zeta_array[0]*1000, zeta_array[2]*1000, color=color, alpha=0.3)
        plt.scatter(np.log10(c_data_fig1[salt]), zeta_data_fig1[salt], label=salt, color=color)
        fontsize = 15
        plt.title('Scales fig 1 - ' + salt, fontsize=fontsize)
        plt.xlabel('log$_{10}(c)$ [mol/L]', fontsize=fontsize)
        plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
        plt.text(-3.5, -120, 'pH = 5.8$\pm$0.3', fontsize=fontsize)
        plt.legend(title='p$K_M$ = %.1f\np$K_a$ = %.1f\n$C_1$ = %i $\mu$F/cm$^2$\n$C_2$ = %i $\mu$F/cm$^2$' % (np.log10(K_ads), pKa, C1*100, C2*100))
        i += 1

# Scales fig 2
plot_scales_fig2 = True
if plot_scales_fig2:
    plt.figure('Scales fig 2')
    plt.scatter(pH_data_fig2, zeta_data_fig2, label='exp')
    pH_list = np.linspace(4, 10, 50)
    zeta_list = np.zeros(len(pH_list))
    i = 0
    for pH in pH_list:
        c_list1 = np.array([1e-3+10**-pH, 1e-3, 10**-pH])
        K_list = np.array([K_K, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_K, C_2=C2_K)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(pH_list, zeta_list*1000, label='model')
    fontsize = 15
    plt.title('Scales fig 2', fontsize=fontsize)
    plt.xlabel('pH', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.legend(title='p$K_M$ = %.1f\np$K_a$ = %.1f\n$C_1$ = %i $\mu$F/cm$^2$\n$C_2$ = %i $\mu$F/cm$^2$' % (np.log10(K_K), pKa, C1_K*100, C2_K*100))

# Osman fig 1
plot_osman_fig1 = True
if plot_osman_fig1:
    pH = 5.8
    cCl = 0.67976
    plt.figure('Osman fig 1')
    plt.scatter(frac_A_list_data_fig1, frac_ads_list_data_fig1)
    frac_A_list = np.linspace(0.1, 0.99, 5)
    frac_ads_list = np.zeros(len(frac_A_list))
    i = 0
    for frac_A in frac_A_list:
        c1 = frac_A*cCl
        c2 = (1-frac_A)*cCl
        c_list1 = np.array([cCl+10**-pH, c1, c2, 10**-pH])
        K_list = np.array([K_Na, K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1, 1])
        v_list = np.array([False, True, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Na, C_2=C2_Na)
        sol.bound_diffuse()
        ads_tot = np.sum(sol.SM_list) + np.sum(sol.bound_diffuse_list[1:-1])
        frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot
        i += 1
    plt.plot(frac_A_list, frac_ads_list)

# Osman, fig 2
plot_osman_fig2 = False
if plot_osman_fig2:
    pH = 5.8
    cCl = 0.6435
    plt.figure('Osman fig 2')
    plt.scatter(frac_A_list_data_fig2, frac_ads_list_data_fig2)
    frac_A_list = np.linspace(0.1, 0.99, 5)
    frac_ads_list = np.zeros(len(frac_A_list))
    i = 0
    for frac_A in frac_A_list:
        c1 = frac_A*cCl
        c2 = (1-frac_A)*cCl
        c_list1 = np.array([cCl+10**-pH, c1, c2, 10**-pH])
        K_list = np.array([K_K, K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1, 1])
        v_list = np.array([False, True, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Na, C_2=C2_Na)
        sol.bound_diffuse()
        ads_tot = np.sum(sol.SM_list) + np.sum(sol.bound_diffuse_list[1:-1])
        frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot
        i += 1
    plt.plot(frac_A_list, frac_ads_list)


######################
# POTENTIAL PROFILES #
######################

# LiCl, vary concentration
plot_LiCl_conc = False
if plot_LiCl_conc:
    pH = 5.8
    plt.figure('LiCl potential profile vary concentration')
    c_list = np.array([1e-1, 1e-2, 1e-3])
    sigma_d_LiCl_conc = np.zeros(len(c_list))
    i = 0
    for c in c_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_PB()
        sigma_d_LiCl_conc[i] = sol.sigma_d
        x = np.concatenate((np.array([-(R_Li+2*R_Li_hyd), -2*R_Li_hyd]), sol.x))
        psi = np.concatenate((np.array([sol.psi_0, sol.psi_beta]), sol.psi))
        plt.plot(x*1e9, psi*1000, label='10$^{%i}$' % np.log10(c))
        i += 1
    plt.title('$\psi$ vs x for LiCl at various conc')
    plt.xlabel('x [nm]')
    plt.ylabel('$\psi$ [mV]')
    plt.legend(title='Conc [M]')
    plt.text(20, -200, 'p$K_M = %.1f$\np$K_a$ = %.1f\npH = %.1f' % (np.log10(K_Li), pKa, pH))

# LiCl, vary K
plot_LiCl_K = False
if plot_LiCl_K:
    c = 1e-1
    pH = 5.8
    plt.figure('LiCl potential profile vary K')
    K_list = np.array([1e2, 1e1, 1e0, 1e-1])
    sigma_d_LiCl_K = np.zeros(len(K_list))
    i = 0
    for K in K_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_PB()
        sigma_d_LiCl_K[i] = sol.sigma_d
        x = np.concatenate((np.array([-(R_Li+2*R_Li_hyd), -2*R_Li_hyd]), sol.x))
        psi = np.concatenate((np.array([sol.psi_0, sol.psi_beta]), sol.psi))
        plt.plot(x*1e9, psi*1000, label='%.1f' % np.log10(K))
        i += 1
    plt.title('$\psi$ vs x for LiCl at various K')
    plt.xlabel('x [nm]')
    plt.ylabel('$\psi$ [mV]')
    plt.legend(title='p$K_M$')
    plt.text(3, -200, 'c = %.1f M\np$K_a$ = %.1f\npH = %.1f' % (c, pKa, pH))


plt.show()
