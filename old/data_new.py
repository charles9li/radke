from old.radii import *
from scipy.constants import epsilon_0
from solver import Force
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

# pKa
pka = 5.3

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

# plot
data_all = False
data_israelachvili = True

##################
# ISRAELACHIVILI #
##################

if data_israelachvili or data_all:
    pH = 5.8
    pKa = 5.3
    c_list1 = ['1e-4', '1e-3', '1e-2', '1e-1']
    c_list1 = ['1e-1']
    for c in c_list1:
        filepath = 'data/israelachvili_fig3_data_' + str(c) + '.csv'
        df_exp = pd.read_csv(filepath_or_buffer=filepath,
                             names=('D', 'F/R'))
        D_list = df_exp.get('D')
        D_start = D_list[0]*1e-9

        c_list = [float(c) + 10**-pH, float(c), 10**-pH]
        K_list = [K_K, 10**pKa]
        z_list = [-1, 1, 1]
        v_list = [False, True, True]
        force = Force(c_list, K_list, z_list, v_list, C1=C1_K, C2=C2_K,
                      pH_effect=False, cation="K")
        force.compute_F(D_start)
        df_model = pd.DataFrame({'D': force.D_list*1e9,
                                 'F/R': force.F_list*1e6})
        df_model.to_csv(path_or_buf='data_figures/israelachvili_'+str(c)+'_model.csv',
                        index=False)
