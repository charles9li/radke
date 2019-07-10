from solver import Force
import pandas as pd

# Parameters
c = 5e-4
pH = 8.9
pKa = 5.3

# Compute force
cH = 10 ** -pH
c_list = [c + cH, c, cH]
K_list = [None, 10 ** pKa]
z_list = [-1, 1, 1]
v_list = [False, True, True]
force = Force(c_list, K_list, z_list, v_list,
              pH_effect=False, cation='Na')
force.compute_F(0.1e-9)
df = pd.DataFrame({'D': force.D_list * 1e9,
                   'F/R': force.F_list * 1e6})
filepath = "data/pashley_fig13_model.csv"
df.to_csv(path_or_buf=filepath, index=False)
