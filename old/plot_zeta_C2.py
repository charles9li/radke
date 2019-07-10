import itertools
from old.solver_old import *

# Parameters
c = 1e-3
pH = 5.8
K_ads = 100
pKa = 5.3
C1 = 10
C2_list = np.linspace(.1, .5, 5)

marker = itertools.cycle(('o', 'v', 's', '8', 'p', '^', '<', '>'))
line = itertools.cycle(("-", "--", ":", "-."))

psi_beta_list = np.zeros(len(C2_list))
zeta_list = np.zeros(len(C2_list))
for C2 in C2_list:
    c_H = 10**-pH
    c_Cl = c + c_H
    c_list = [c_Cl, c, c_H]
    K_list = np.array([K_ads, 10**pKa])
    z_list = np.array([-1, 1, 1])
    v_list = np.array([False, True, True])
    sol = Solution_1plate(c_list, K_list, z_list, v_list, pH_effect=False, C_1=C1, C_2=C2)
    sol.solver_sigma()
    zeta_list[i] =
