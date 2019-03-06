import matplotlib.pyplot as plt
from data_scales import *
from data_osman import *
from radii import *
from solver import *
from line_marker import *
from scipy.constants import epsilon_0


##############
# PARAMETERS #
##############

# K_ads
K_Li = 10**0.1
K_Na = 10**0.1
K_K = 10**2.5
K_Cs = 10**3
K_ads_list = [K_Li, K_Na, K_K, K_Cs]

# pH
pH = 5.8

# pKa
pKa = 5.6

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
plt.figure('Scales fig 1')
i = 0
for salt in salts:
    K_ads = K_ads_list[i]
    C1 = C1_list[i]
    C2 = C2_list[i]
    c_list = 10**np.linspace(-5, -1, 100)
    zeta_list = np.zeros(len(c_list))
    j = 0
    for c in c_list:
        c_list1 = np.array([c+10**-pH, c, 10**-pH])
        K_list = np.array([K_ads, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list1, K_list, z_list, v_list,
                              pH_effect=False, C_1=C1, C_2=C2)
        sol.solver_sigma()
        zeta_list[j] = sol.psi_d
        j += 1
    i += 1
    plt.scatter(np.log10(c_data_fig1[salt]), zeta_data_fig1[salt],
                marker=next(markers), label=salt)
    plt.plot(np.log10(c_list), zeta_list*1000,
             linestyle=next(lines))
font_size = 15
plt.rc('font', size=font_size)
plt.rc('axes', labelsize=font_size)
plt.rc('xtick', labelsize=font_size)
plt.rc('ytick', labelsize=font_size)
plt.title('Scales fig 1')
plt.xlabel('log$_{10}(c)$ [mol/L]', fontsize=font_size)
plt.ylabel('Zeta potential [mV]', fontsize=font_size)
plt.legend()

plt.show()
