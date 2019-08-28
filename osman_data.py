from solver import Solution1Plate
import numpy as np
import pandas as pd

# Parameters
V = 0.025
pH = 5.8
pKa = 5.3

# Lists
fig_index_list = [1, 2, 4]
cA0_list = [0.64, 0.62, 0.67]
A_list = ["Na", "K", "Cs"]

# Process data
for i in range(len(fig_index_list)):

    # Extract data from csv file
    fig_index = fig_index_list[i]
    filepath = "data/osman fig{0} data.csv".format(fig_index)
    df = pd.read_csv(filepath_or_buffer=filepath,
                     names=("A0", "GammaA", "cLi"))

    # Store data into variables
    A0_list = np.array(df.get("A0"))
    GammaA_list = np.array(df.get("GammaA"))
    cLi_list = np.array(df.get("cLi"))

    # Compute mass of mica
    cA0 = cA0_list[i]
    m_list = 1000*cA0*V/A0_list

    # Compute cA
    cA_list = (cA0*V-GammaA_list*m_list/1000)/V

    # Prune lists
    j = 0
    while cA_list[j] < 0:
        j += 1
    A0_list = A0_list[j:]
    GammaA_list = GammaA_list[j:]
    cLi_list = cLi_list[j:]
    m_list = m_list[j:]

    # Change cLi to mM
    cLi_list = cLi_list*m_list/(1000*V)

    # Compute exchange isotherms
    SA_list = np.zeros(len(cA_list))
    SLi_list = np.zeros(len(cA_list))
    SH_list = np.zeros(len(cA_list))
    for j in range(len(cA_list)):
        cA = cA_list[j]*1e-3
        cLi = cLi_list[j]*1e-3
        cH = 10**-pH
        cCl = cA + cLi + cH
        c_list = [cCl, cA, cLi, cH]
        K_list = [None, None, 10**5.3]
        z_list = [-1, 1, 1, 1]
        v_list = [False, True, True, True]
        sol = Solution1Plate(c_list, K_list, z_list, v_list, pH_effect=False,
                             cation=[A_list[i], "Li"])
        sol.solve_equations()
        SA_list[j] = sol.SM_list[0]
        SLi_list[j] = sol.SM_list[1]
        SH_list[j] = sol.SM_list[2]
