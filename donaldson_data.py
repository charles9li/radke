from solver import Force
import pandas as pd

# Parameters
pKa = 5.3

# Lists
pH_list = [3, 10]
c_list = [1, 5, 10]

# Compute force
for pH in pH_list:
    for c in c_list:
        cH = 10 ** -pH
        c_list1 = [c * 1e-3 + cH, c * 1e-3, cH]
        K_list = [None, 10 ** pKa]
        z_list = [-1, 1, 1]
        v_list = [False, True, True]
        force = Force(c_list1, K_list, z_list, v_list,
                      pH_effect=False, cation='Na')
        force.compute_F(0.1e-9)

        filepath = 'data/donaldson_pH' + str(pH) + '_' + str(c) + 'mM_model.csv'
        df = pd.DataFrame({'D': force.D_list * 1e9,
                           'F/R': force.F_list * 1e6})
        df.to_csv(path_or_buf=filepath,
                  index=False)
