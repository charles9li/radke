from solver import Force
import pandas as pd

# Parameters
pH = 5.8
del_pH = 0.3
pKa = 5.3

# Lists
c_list = ['1e-4', '1e-3', '1e-2', '1e-1']

# Produce data
for c in c_list:
    c_list_pH_low = [float(c) + 10 ** -(pH - del_pH), float(c), 10 ** -(pH - del_pH)]
    c_list_pH_mid = [float(c) + 10 ** -pH, float(c), 10 ** -pH]
    c_list_pH_high = [float(c) + 10 ** -(pH + del_pH), float(c), 10 ** -(pH + del_pH)]
    K_list = [None, 10 ** pKa]
    z_list = [-1, 1, 1]
    v_list = [False, True, True]

    force_pH_low = Force(c_list_pH_low, K_list, z_list, v_list,
                         pH_effect=False, cation='K')
    force_pH_low.compute_F(0.1e-9)

    force_pH_mid = Force(c_list_pH_mid, K_list, z_list, v_list,
                         pH_effect=False, cation='K')
    force_pH_mid.compute_F(0.1e-9)

    force_pH_high = Force(c_list_pH_high, K_list, z_list, v_list,
                          pH_effect=False, cation='K')
    force_pH_high.compute_F(0.1e-9)

    length = min([len(force_pH_low.D_list), len(force_pH_mid.D_list), len(force_pH_high.D_list)])
    df = pd.DataFrame({'D': force_pH_low.D_list[:length] * 1e9,
                       'F/R_pH_low': force_pH_low.F_list[:length] * 1e6,
                       'F/R_pH_mid': force_pH_mid.F_list[:length] * 1e6,
                       'F/R_pH_high': force_pH_high.F_list[:length] * 1e6})
    filepath = 'data/israelachvili_' + c + '_model.csv'
    df.to_csv(path_or_buf=filepath,
              index=False)
