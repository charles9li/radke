from radii import *
from scipy.constants import epsilon_0, k
from solver import Solution_1plate
from solvers import *
import numpy as np
import pandas as pd


##############
# PARAMETERS #
##############

# K_ads
K_Li = 10**0.3
K_Na = K_Li*2.5
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

# compute all data
data_all = False


##########
# SCALES #
##########

# Figure 1
data_scales_fig1 = True
if data_scales_fig1 or data_all:
    salts = ['LiCl', 'NaCl', 'KCl', 'CsCl']
    i = 0
    for salt in salts:
        K_ads = K_ads_list[i]
        C1 = C1_list[i]
        C2 = C2_list[i]
        c_list = np.logspace(-5, -1, 45)
        zeta_array = np.zeros((3, len(c_list)))
        j = 0
        for pH in np.array([5.5, 5.8, 6.1]):
            k = 0
            for c in c_list:
                c_list1 = np.array([c+10**-pH, c, 10**-pH])
                K_list = np.array([K_ads, 10**pKa])
                z_list = np.array([-1, 1, 1])
                v_list = np.array([False, True, True])
                sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                                      pH=pH, pH_effect=False, C_1=C1, C_2=C2)
                sol.solver_sigma()
                zeta_array[j, k] = sol.psi_d*1000
                k += 1
            j += 1
        df = pd.DataFrame({'c': c_list,
                           'zeta_pH5.5': zeta_array[0],
                           'zeta_pH5.8': zeta_array[1],
                           'zeta_pH6.1': zeta_array[2]})
        df.to_csv(path_or_buf='data_figures/scales_fig1_model_'+salt+'.csv',
                  index=False)
        i += 1

# Figure 2
data_scales_fig2 = False
if data_scales_fig2 or data_all:
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
    df = pd.DataFrame({'pH': pH_list,
                       'zeta': zeta_list*1000})
    df.to_csv(path_or_buf='data_figures/scales_fig2_model.csv',
              index=False)


#########
# OSMAN #
#########

data_osman = False
if data_osman or data_all:
    cCl_list = [0.67976, 0.6435, 0.66163]
    Gamma_Li0_list = [459.74, 468.8737, 608.0727]
    CA0_list = [0.64, 0.62, 0.67]
    K_list_osman = [K_Na, K_K, K_Cs]
    C1_list_osman = [C1_Na, C1_K, C1_Cs]
    C2_list_osman = [C2_Na, C2_K, C2_Cs]
    i = 0
    for fig_index in [1, 2, 4]:
        cCl = cCl_list[i]
        Gamma_Li0 = Gamma_Li0_list[i]
        CA0 = CA0_list[i]
        df = pd.read_csv(filepath_or_buffer='data/osman fig'+str(fig_index)+' data.csv',
                         names=('A0', 'Gamma_A', 'cLi'))
        Gamma_Li_list = Gamma_Li0*np.ones(len(df.get('cLi'))) - df.get('cLi')
        cA_list = CA0*(1-np.divide(df.get('Gamma_A'), df.get('A0')))
        df_exp = pd.DataFrame({'frac_A': cA_list/cCl,
                               'frac_ads': df.get('Gamma_A')/(df.get('Gamma_A')+Gamma_Li_list)})
        df_exp.to_csv(path_or_buf='data_figures/osman_fig'+str(fig_index)+'_exp.csv',
                      index=False)

        frac_A_list = np.linspace(0.0001, 0.1, 20)
        frac_A_list = np.append(frac_A_list, np.linspace(0.11, 0.99, 15))
        frac_ads_list = np.zeros(len(frac_A_list))
        j = 0
        for frac_A in frac_A_list:
            K_ads = K_list_osman[i]
            C1 = C1_list_osman[i]
            C2 = C2_list_osman[i]
            c1 = frac_A*cCl
            c2 = (1-frac_A)*cCl
            c_list = np.array([cCl+10**-pH, c1, c2, 10**-pH])
            K_list = np.array([K_ads, K_Li, 10**pKa])
            z_list = np.array([-1, 1, 1, 1])
            v_list = np.array([False, True, True, True])
            sol = Solution_1plate(c_list, K_list, z_list, v_list,
                                  pH_effect=False, C_1=C1, C_2=C2)
            # sol.bound_diffuse()
            # ads_tot = np.sum(sol.SM_list[0:2]) + np.sum(sol.bound_diffuse_list[1:-1])
            # frac_ads_list[j] = (sol.SM_list[0] + sol.bound_diffuse_list[1])/ads_tot

            sol.solver_sigma()
            frac_ads_list[j] = sol.SM_list[0]/(np.sum(sol.SM_list[0:2]))
            j += 1

        df_model = pd.DataFrame({'frac_A': frac_A_list,
                                 'frac_ads': frac_ads_list})
        df_model.to_csv(path_or_buf='data_figures/osman_fig'+str(fig_index)+'_model.csv',
                        index=False)

        i += 1


#################
# ISRAELACHVILI #
#################

def compute_kappa(c, pH=5.8, T=298, epsilon=79):
    rho = 1000*N_A*(c+10**-pH)
    return np.sqrt(2*rho*e**2/(epsilon*epsilon_0*k*T))

data_israelachvili = False
if data_israelachvili or data_all:
    pH = 5.8
    c_list = ['1e-4', '1e-3', '1e-2', '1e-1']
    for c in c_list:
        df_exp = pd.read_csv(filepath_or_buffer='data/israelachvili_fig3_data_'+str(c)+'.csv',
                         names=('D', 'F/R'))
        D_list = df_exp.get('D')
        D_list = np.linspace(D_list[0], D_list[len(D_list)-1], 50)*1e-9
        FR_list = np.zeros(len(D_list))
        i = 0
        for D in D_list:
            sol = solver_2plate_2cation_ads(float(c), 10**-pH,
                                            K_K, 10**pKa, D,
                                            C_1=C1_K, C_2=C2_K)
            rho_c_interp = interp1d(sol.x, sol.cation_1, 'cubic')
            rho_H_interp = interp1d(sol.x, sol.cation_2, 'cubic')
            rho_Cl_interp = interp1d(sol.x, sol.anion, 'cubic')
            rho_tot_m = rho_c_interp(D/2) + rho_H_interp(D/2) + rho_Cl_interp(D/2)
            rho_tot_bulk = 2*N_A*1000*(float(c) + 10**-pH)

            P = k*298*(rho_tot_m - rho_tot_bulk)
            W = P/compute_kappa(float(c))
            FR = 2*np.pi*W
            W_H = -2.2e-20/(12*np.pi*D**2)
            FR_H = 2*np.pi*W_H
            FR_list[i] = FR + FR_H
            i += 1
        df_model = pd.DataFrame({'D': D_list*1e9,
                                 'F/R': FR_list*1e6})
        df_model.to_csv(path_or_buf='data_figures/israelachvili_'+str(c)+'_model.csv',
                        index=False)


#########
# SIGMA #
#########

data_sigma_beta_conc = False
if data_sigma_beta_conc:
    pH = 5.8
    c_list = np.logspace(-5, -1, 50)
    sigma_beta_frac_list = np.zeros(len(c_list))
    i = 0
    for c in c_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_sigma()
        sigma_beta_frac_list[i] = np.sum(sol.sigma_beta_list)/-sol.sigma_0
        i += 1
    df = pd.DataFrame({'c': c_list,
                       'beta_frac': sigma_beta_frac_list})
    df.to_csv(path_or_buf='data_figures/sigma_beta_conc.csv',
              index=False)

#
data_sigma_beta_H = True
if data_sigma_beta_H:
    pH = 5.8
    c_list = np.logspace(-5, -1, 50)
    sigma_beta_frac_M_list = np.zeros(len(c_list))
    sigma_beta_frac_H_list = np.zeros(len(c_list))
    i = 0
    for c in c_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K_Li, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1_Li, C_2=C2_Li)
        sol.solver_sigma()
        sigma_beta_frac_M_list[i] = sol.sigma_beta_list[0]/-sol.sigma_0
        sigma_beta_frac_H_list[i] = sol.sigma_beta_list[-1]/-sol.sigma_0
        i += 1
    df = pd.DataFrame({'c': c_list,
                       'beta_frac_M': sigma_beta_frac_M_list,
                       'beta_frac_H': sigma_beta_frac_H_list})
    df.to_csv(path_or_buf='data_figures/sigma_beta_conc_H.csv',
              index=False)