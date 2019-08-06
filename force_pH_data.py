from solver import Force
import pandas as pd

# Parameters
c = 1e-3
cation = 'K'
pKa = 5.3

# Lists
pH_list = [4, 5, 6, 8, 10]

# Compute force curves
for pH in pH_list:
    cH = 10 ** -pH
    c_list = [c + cH, c, cH]
    K_list = [None, 10 ** pKa]
    z_list = [-1, 1, 1]
    v_list = [False, True, True]
    force = Force(c_list, K_list, z_list, v_list, pH_effect=False, cation=cation)
    force.compute_F(0.1e-9)

    df = pd.DataFrame({'D': force.D_list * 1e9,
                       'F/R': force.F_list * 1e6})
    filepath = 'data/force_' + cation + 'Cl_c' + str(c) + '_pH' + str(pH) + '_model.csv'
    df.to_csv(path_or_buf=filepath, index=False)
