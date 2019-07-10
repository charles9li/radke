from solver import Force
import pandas as pd

# Parameters
pKa = 5.3

# Lists
fig_index_list = [1, 2, 3, 4]
cation_list = ['Li', 'Na', 'K', 'Cs']
c_list_list = [['6e-2', '1e-2', '1e-3', '1e-4'],
               ['1e-2', '1e-3', '1e-4', '4e-5'],
               ['1e-3', '3e-4', '4e-5'],
               ['1e-3', '1e-4', '4e-5']]

pH_list = [5.4, 5.7, 5.7, 5.8]

# Compute force
for i in range(len(fig_index_list)):
    fig_index = str(fig_index_list[i])
    cation = cation_list[i]
    c_list = c_list_list[i]
    pH = pH_list[i]
    for c in c_list:
        c_list_pH_low = [float(c) + 10 ** -(pH-0.3), float(c), 10 ** -(pH-0.3)]
        c_list_pH_mid = [float(c) + 10 ** -pH, float(c), 10 ** -pH]
        c_list_pH_high = [float(c) + 10 ** -(pH+0.3), float(c), 10 ** -(pH+0.3)]
        K_list = [None, 10**5.3]
        z_list = [-1, 1, 1]
        v_list = [False, True, True]
        force_pH_low = Force(c_list_pH_low, K_list, z_list, v_list,
                             pH_effect=False, cation=cation)
        force_pH_low.compute_F(0.1e-9)
        force_pH_mid = Force(c_list_pH_mid, K_list, z_list, v_list,
                             pH_effect=False, cation=cation)
        force_pH_mid.compute_F(0.1e-9)
        force_pH_high = Force(c_list_pH_high, K_list, z_list, v_list,
                              pH_effect=False, cation=cation)
        force_pH_high.compute_F(0.1e-9)
        length = min([len(force_pH_low.D_list), len(force_pH_mid.D_list), len(force_pH_high.D_list)])
        df_model = pd.DataFrame({'D': force_pH_mid.D_list[:length] * 1e9,
                                 'F/R_pH_low': force_pH_low.F_list[:length] * 1e6,
                                 'F/R_pH_mid': force_pH_mid.F_list[:length] * 1e6,
                                 'F/R_pH_high': force_pH_high.F_list[:length] * 1e6})
        filepath = "data/pashley_fig" + fig_index + "_model_" + c + ".csv"
        df_model.to_csv(path_or_buf=filepath, index=False)
