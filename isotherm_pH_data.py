from solver import Solution1Plate
import numpy as np
import pandas as pd

# Parameters
pH = 5.8
pKa = 5.3

# Lists
cation_list = ["Na", "K", "Cs"]
concentration_list = ["1e-3", "0.6"]

# Compute isotherms
sol_prev = None
for concentration in concentration_list:
    c = float(concentration)
    for cation in cation_list:
        fraction_list = np.linspace(0.001, 0.99, 50)
        gamma_M_list = np.zeros(len(fraction_list))
        gamma_Li_list = np.zeros(len(fraction_list))
        gamma_H_list = np.zeros(len(fraction_list))
        for i in range(len(fraction_list)):
            fraction = fraction_list[i]
            c_list = [c, fraction * c, (1 - fraction) * c]
            K_list = [None, None]
            z_list = [-1, 1, 1]
            v_list = [False, True, True]
            sol = Solution1Plate(c_list, K_list, z_list, v_list,
                                 pH_effect=True, pH=pH, pKa=pKa,
                                 cation=[cation, 'Li'])
            sol.solve_equations(sol_prev)

            gamma_M_list[i] = sol.SM_list[0]/sol.L
            gamma_Li_list[i] = sol.SM_list[1]/sol.L
            gamma_H_list[i] = sol.SH/sol.L
        df = pd.DataFrame({"fraction": fraction_list,
                           "gamma_M": gamma_M_list,
                           "gamma_Li": gamma_Li_list,
                           "gamma_H": gamma_H_list})
        filepath = "data_figures/isotherm_pH_" + cation + "_" + concentration + "M_model.csv"
        df.to_csv(path_or_buf=filepath, index=False)
