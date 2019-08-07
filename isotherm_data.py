from solver import Solution1Plate
import numpy as np
import pandas as pd

# Parameters
c = 0.6
pH = 5.8
pKa = 5.3

# Lists
cation_list = ["K", "Cs"]

# Compute isotherms
sol_prev = None
for cation in cation_list:
    fraction_list = np.linspace(0.001, 0.99, 50)
    gamma_M_list = np.zeros(len(fraction_list))
    gamma_Li_list = np.zeros(len(fraction_list))
    gamma_H_list = np.zeros(len(fraction_list))
    for i in range(len(fraction_list)):
        fraction = fraction_list[i]
        c_list = [c + 10 ** -pH, fraction * c, (1 - fraction) * c, 10 ** -pH]
        K_list = [None, None, 10 ** pKa]
        z_list = [-1, 1, 1, 1]
        v_list = [False, True, True, True]
        sol = Solution1Plate(c_list, K_list, z_list, v_list,
                             pH_effect=False,
                             cation=[cation, 'Li'])
        sol.solve_equations(sol_prev)

        gamma_M_list[i] = sol.SM_list[0]/sol.L
        gamma_Li_list[i] = sol.SM_list[1]/sol.L
        gamma_H_list[i] = sol.SM_list[2]/sol.L
    df = pd.DataFrame({"fraction": fraction_list,
                       "gamma_M": gamma_M_list,
                       "gamma_Li": gamma_Li_list,
                       "gamma_H": gamma_H_list})
    filepath = "data_figures/isotherm_" + cation + "_model.csv"
    df.to_csv(path_or_buf=filepath, index=False)
