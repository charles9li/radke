from scipy.constants import e, epsilon_0, k, N_A
from solver import Solution2Plate
import numpy as np
import pandas as pd

# Parameters
c = 1e-3
pH = 5.8
pKa = 5.3
eps = 79 * epsilon_0
T = 298

# Lists
cation_list = ['Li', 'Na', 'K', 'Cs']

# Initialize prev guess
sol_prev = None

# Compute data
i = 0
for cation in cation_list:
    D_list = np.linspace(0.1, 50, 50)*1e-9
    sigma_d_list = np.zeros(len(D_list))
    sigma_d_grahame_list = np.zeros(len(D_list))
    j = 0
    for D in D_list:
        cH = 10 ** -pH
        c_list = [c + cH, c, cH]
        K_list = [None, 10 ** pKa]
        z_list = [-1, 1, 1]
        v_list = [False, True, True]
        sol = Solution2Plate(c_list, K_list, z_list, v_list, D,
                             pH_effect=False, cation=cation)
        sol.solve_equations(sol_prev)
        sol_prev = sol
        sigma_d_list[j] = sol.sigma_d
        psi_d = sol.psi_profile(0)
        sigma_d_grahame_list[j] = -np.sqrt(8*eps*k*T*N_A*1000*c)*np.sinh(e*psi_d/(2*k*T))
        j += 1
    df = pd.DataFrame({'D': D_list*1e9,
                       'sigma_d': sigma_d_list,
                       'sigma_d_grahame': sigma_d_grahame_list})
    df.to_csv(path_or_buf='data/grahame_' + cation + 'Cl.csv',
              index=False)
    i += 1
