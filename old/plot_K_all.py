from old.data_claesson import *
from old.data_israelachvili import *
from old.data_osman import *
from old.data_pashley import *
from old.data_scales import *
from line_marker import *
from old.radii import *
from old.solver_old import *
from old.solvers import *


##############
# PARAMETERS #
##############

# K_ads
K_Li = 10**0.3
K_Na = 10**0.7
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


###################
# PLOT PERMISSION #
###################

plot_scales_fig1 = False
plot_scales_fig2 = False
plot_osman_fig1 = True
plot_osman_fig2 = False
plot_osman_fig4 = False
plot_claesson_fig2 = False
plot_sigma_K = False
plot_LiCl_conc = False
plot_LiCl_K = False
plot_pashley_fig1 = False
plot_pashley_fig2 = False
plot_pashley_fig3 = False
plot_pashley_fig3_inset = False
plot_pashley_fig4 = False
plot_israelachvili_fig3 = False


#########
# PLOTS #
#########

# Scales fig 1
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
if plot_osman_fig1:
    pH = 5.8
    cCl = 0.67976
    plt.figure('Osman fig 1')
    plt.scatter(frac_A_list_data_fig1, frac_ads_list_data_fig1)
    frac_A_list = np.linspace(0.1, 0.99, 50)
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
        ads_tot = np.sum(sol.SM_list[0:2]) + np.sum(sol.bound_diffuse_list[1:-1])
        frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot

        # sol.solver_sigma()
        # frac_ads_list[i] = sol.SM_list[0]/(np.sum(sol.SM_list[0:2]))

        i += 1
    plt.plot(frac_A_list, frac_ads_list)
    fontsize = 15
    plt.title('Osman Fig 1 - Na-Li', fontsize=fontsize)
    plt.xlabel('[Na+]/[Cl-]', fontsize=fontsize)
    plt.ylabel('$\Gamma_M$/CEC', fontsize=fontsize)

# Osman, fig 2
if plot_osman_fig2:
    pH = 5.8
    cCl = 0.6435
    plt.figure('Osman fig 2')
    plt.scatter(frac_A_list_data_fig2, frac_ads_list_data_fig2)
    frac_A_list = np.linspace(0.001, 0.99, 100)
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
                              pH_effect=False, C_1=C1_K, C_2=C2_K)
        sol.bound_diffuse()
        ads_tot = np.sum(sol.SM_list) + np.sum(sol.bound_diffuse_list[1:-1])
        frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot

        # sol.solver_sigma()
        # frac_ads_list[i] = sol.SM_list[0]/(np.sum(sol.SM_list[0:2]))

        i += 1
    plt.plot(frac_A_list, frac_ads_list)
    fontsize = 15
    plt.title('Osman Fig 2 - K-Li', fontsize=fontsize)
    plt.xlabel('[K+]/[Cl-]', fontsize=fontsize)
    plt.ylabel('$\Gamma_M$/CEC', fontsize=fontsize)

# Osman, fig 4
if plot_osman_fig4:
    pH = 5.8
    cCl = 0.66163
    plt.figure('Osman fig 4')
    plt.scatter(frac_A_list_data_fig4, frac_ads_list_data_fig4)
    frac_A_list = np.linspace(0.001, 0.99, 100)
    frac_ads_list = np.zeros(len(frac_A_list))
    i = 0
    for frac_A in frac_A_list:
        c1 = frac_A*cCl
        c2 = (1-frac_A)*cCl
        c_list1 = np.array([cCl+10**-pH, c1, c2, 10**-pH])
        K_list = np.array([K_Cs, K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1, 1])
        v_list = np.array([False, True, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Na, C_2=C2_Na)
        sol.bound_diffuse()
        ads_tot = np.sum(sol.SM_list) + np.sum(sol.bound_diffuse_list[1:-1])
        frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot

        # sol.solver_sigma()
        # frac_ads_list[i] = sol.SM_list[0]/(np.sum(sol.SM_list[0:2]))

        i += 1
    plt.plot(frac_A_list, frac_ads_list)
    fontsize = 15
    plt.title('Osman Fig 4 - Cs-Li', fontsize=fontsize)
    plt.xlabel('[Cs+]/[Cl-]', fontsize=fontsize)
    plt.ylabel('$\Gamma_M$/CEC', fontsize=fontsize)

# Claesson, fig 2
if plot_claesson_fig2:
    pH = 5.8
    plt.figure('Claesson fig 2')
    plt.scatter(
        claesson_fig2_data.wash[0], claesson_fig2_data.wash[1],
        label='wash')
    plt.scatter(
        claesson_fig2_data.press[0], claesson_fig2_data.press[1],
        label='press')
    M_list = np.linspace(1, 6, 100)
    S_list = np.zeros(len(M_list))
    i = 0
    for M in M_list:
        cH = 10**-pH
        c = M*10*cH
        cCl = c + cH
        c_list = np.array([cCl, c, cH])
        K_list = np.array([K_Na, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list, K_list, z_list, v_list,
                              pH=pH, pKa=pKa, pH_effect=False,
                              C_1=C1_Na, C_2=C2_Na)
        sol.solver_sigma()
        S_list[i] = sol.SM_list[0]/100
        i += 1
    plt.semilogy(M_list, S_list, label='model')
    plt.legend()
    plt.yscale('log')
    plt.xlim((7, 0))


#########
# SIGMA #
#########

if plot_sigma_K:
    pH = 5.8
    pKa = 3
    plt.figure('sigma vary K c = 1e-3')
    c_list = np.array([1e-1, 1e-2, 1e-3])
    for c in c_list:
        K_list = 10**np.linspace(-3, 3, 20)
        sigma_d_list = np.zeros(len(K_list))
        i = 0
        for K in K_list:
            c_list1 = np.array([c+10**-pH, c, 10**-pH])
            K_list1 = np.array([K, 10**pKa])
            z_list = np.array([-1, 1, 1])
            v_list = np.array([False, True, True])
            sol = Solution_1plate(c_list1, K_list1, z_list, v_list,
                                  pH_effect=False, C_1=C1_Li, C_2=C2_Li)
            sol.solver_sigma()
            sigma_d_list[i] = sol.sigma_d
            i += 1
        plt.plot(np.log10(K_list), sigma_d_list, label='10$^{%i}$' % np.log10(c))
    fontsize = 15
    plt.title('$\sigma_d$ vs $K_M$ for various conc', fontsize=fontsize)
    plt.xlabel('p$K_M$', fontsize=fontsize)
    plt.ylabel('$\sigma_d$ [C/m$^2$]', fontsize=fontsize)
    plt.text(-2, 0.003, 'p$K_a$ = %.1f\npH = %.1f\n$C_1$ = %i $\mu$F/cm$^2$\n$C_2$ = %i $\mu$F/cm$^2$' % (pKa, pH, C1_Li*100, C2_Li*100), fontsize=fontsize)
    plt.legend(title='Conc [M]')


######################
# POTENTIAL PROFILES #
######################

# LiCl, vary concentration
if plot_LiCl_conc:
    pH = 5.8
    plt.figure('LiCl potential profile vary concentration')
    c_list = np.array([1e-1, 1e-2, 1e-3, 1e-4])
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
    plt.plot(np.array([0, 0]), np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd, -2*R_Li_hyd])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    fontsize = 15
    plt.title('$\psi$ vs x for LiCl at various conc', fontsize=fontsize)
    plt.xlabel('x [nm]', fontsize=fontsize)
    plt.ylabel('$\psi$ [mV]', fontsize=fontsize)
    plt.legend(title='Conc [M]')
    plt.text(20, -400, 'p$K_M = %.1f$\np$K_a$ = %.1f\npH = %.1f\n$C_1$ = %i $\mu$F/cm$^2$\n$C_2$ = %i $\mu$F/cm$^2$' % (np.log10(K_Li), pKa, pH, C1_Li*100, C2_Li*100), fontsize=fontsize)

# LiCl, vary K
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
    plt.plot(np.array([0, 0]), np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    plt.plot(np.array([-2*R_Li_hyd, -2*R_Li_hyd])*1e9, np.array([-700, 0]), color='k', linestyle='--', linewidth=1)
    fontsize = 15
    plt.title('$\psi$ vs x for LiCl at various K', fontsize=fontsize)
    plt.xlabel('x [nm]', fontsize=fontsize)
    plt.ylabel('$\psi$ [mV]', fontsize=fontsize)
    plt.legend(title='p$K_M$')
    plt.text(2, -400, 'c = %.1f M\np$K_a$ = %.1f\npH = %.1f\n$C_1$ = %i $\mu$F/cm$^2$\n$C_2$ = %i $\mu$F/cm$^2$' % (c, pKa, pH, C1_Li*100, C2_Li*100), fontsize=fontsize)


################
# FORCE CURVES #
################

def compute_kappa(c, pH=5.8, T=298, epsilon=79):
    rho = 1000*N_A*(c+10**-pH)
    return np.sqrt(2*rho*e**2/(epsilon*epsilon_0*k*T))


# Pashley, fig 1
if plot_pashley_fig1:
    plt.figure('Pashley fig 1')
    plot_indices = np.array([1, 2])
    for plot_index in plot_indices:
        data = pashley_fig1_data[plot_index]
        plt.scatter(data.D, data.FR/1e3)
        D_list = np.linspace(data.D[0], data.D[-1], 50)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_Li, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        plt.semilogy(D_list*1e9, FR_list*1000, label='%.0e' % data.c)
    plt.title('Pashley Fig 1 - LiCl')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.legend(title='Conc [M]')

# Pashley, fig 2
if plot_pashley_fig2:
    plt.figure('Pashley fig 2')
    plot_indices = np.array([1, 2])
    for plot_index in plot_indices:
        data = pashley_fig2_data[plot_index]
        plt.scatter(data.D, data.FR/1e3)
        D_list = np.linspace(data.D[0], data.D[-1], 50)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_Na, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        plt.semilogy(D_list*1e9, FR_list*1000, label='%.0e' % data.c)
    plt.title('Pashley Fig 2 - NaCl')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.legend(title='Conc [M]')

# Pashley, fig 3
if plot_pashley_fig3:
    plt.figure('Pashley fig 3')
    plot_indices = np.array([1, 2])
    for plot_index in plot_indices:
        data = pashley_fig3_data[plot_index]
        plt.scatter(data.D, data.FR/1e3)
        D_list = np.linspace(data.D[0], data.D[-1], 3)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_Na, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        plt.semilogy(D_list*1e9, FR_list*1000, label='%.0e' % data.c)
    plt.title('Pashley Fig 3 - KCl')
    plt.xlabel('$D$ [nm]')
    plt.ylabel('$F/R$ [mN/m]')
    plt.legend(title='Conc [M]')

# Pashley, fig 3, inset
if plot_pashley_fig3_inset:

    # Format figure
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    fig, ax1 = plt.subplots()
    fig.set_size_inches(8, 6)
    left, bottom, width, height = [0.50, 0.55, 0.38, 0.3]
    ax2 = fig.add_axes([left, bottom, width, height])
    left, bottom, width, height = [0.2, 0.195, 0.17, 0.18]
    ax3 = fig.add_axes([left, bottom, width, height])
    # left, bottom, width, height = [0.4, 0.55, 0.48, 0.3]
    # ax2 = fig.add_axes([left, bottom, width, height])
    # left, bottom, width, height = [0.3, 0.6, 0.2, 0.18]
    # ax3 = fig.add_axes([left, bottom, width, height])

    # Plot forces
    plot_indices = np.array([1, 2])

    for plt_index in plot_indices:
        color = next(colors)
        pH = 5.7
        data = pashley_fig3_data[plt_index]
        ax1.scatter(data.D, data.FR/1e3,
                    label=data.c_str,
                    marker=next(markers),
                    color=color)
        D_list = np.linspace(.4, data.D[-1], 100)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_Na, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        ax1.semilogy(D_list*1e9, FR_list*1000,
                     linestyle=next(lines),
                     color=color)
    ax1.set_xlabel('$D$ [nm]')
    ax1.set_ylabel('$F/R$ [mN/m]')
    ax1.legend(
        loc='lower left',
        title='c [M]',
        bbox_to_anchor=(0.1, 0.75))

    # Plot isotherm
    color = next(colors)
    pH = 5.8
    cCl = 0.6435
    ax3.scatter(frac_A_list_data_fig2, frac_ads_list_data_fig2,
                color=color,
                marker=next(markers))
    frac_A_list = np.linspace(0.001, 0.99, 40)
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
                              pH_effect=False, C_1=C1_K, C_2=C2_K)
        # sol.bound_diffuse()
        # ads_tot = np.sum(sol.SM_list) + np.sum(sol.bound_diffuse_list[1:-1])
        # frac_ads_list[i] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot

        sol.solver_sigma()
        frac_ads_list[i] = sol.SM_list[0]/(np.sum(sol.SM_list[0:2]))

        i += 1
    ax3.plot(frac_A_list, frac_ads_list,
             color=color)
    ax3.set_xlabel('[K$^+$] / [Cl$^-$]')
    ax3.set_ylabel('$\Gamma_K$ / CEC')

    # Plot zeta potential
    color = next(colors)
    ax2.scatter(np.log10(c_data_fig1['KCl']), zeta_data_fig1['KCl'],
                color=color,
                marker=next(markers))
    K_ads = K_K
    C1 = C1_K
    C2 = C2_K
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
    ax2.plot(np.log10(c_list), zeta_array[1]*1000,
             color=color)
    ax2.set_xlabel('log($c$) [mol dm$^{-3}$]')
    ax2.set_ylabel('$\zeta$ [mV]')

# Pashley, fig 4
if plot_pashley_fig4:
    plt.figure('Pashley fig 4')
    plot_indices = np.array([2])
    for plot_index in plot_indices:
        data = pashley_fig4_data[plot_index]
        plt.scatter(data.D, data.FR/1e3)
        D_list = np.linspace(data.D[0], data.D[-1], 50)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_Na, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        plt.semilogy(D_list*1e9, FR_list*1000, label='%.0e' % data.c)
    plt.title('Pashley Fig 4 - CsCl')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.legend(title='Conc [M]')

# Israelachvili, fig 3
if plot_israelachvili_fig3:
    plt.figure('Israelachvili fig 3')
    plot_indices = np.array([0, 1, 2, 3])
    for plot_index in plot_indices:
        data = israelachvili_fig3_data[plot_index]
        plt.scatter(data.D, data.FR/1e3)
        D_list = np.linspace(data.D[0], data.D[-1], 50)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(data.c, 10**-pH, K_K, 10**pKa, D,
                                            C_1=C1_Li, C_2=C2_Li)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(data.c + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(data.c)
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1

        plt.semilogy(D_list*1e9, FR_list*1000, label='%.0e' % data.c)
    plt.title('Israelachvili Fig 3 - KNO$_3$')
    plt.xlabel('D [nm]')
    plt.ylabel('F/R [mN/m]')
    plt.legend(title='Conc [M]')


plt.show()
