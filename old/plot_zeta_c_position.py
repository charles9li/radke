import matplotlib.pyplot as plt
import itertools
from solver import *
from scipy.interpolate import interp1d

# Parameters
log10c_list = np.linspace(-5, -1, 10)
K_ads = 10**3.151010098187411
pH = 5.8
pKa = 5.3
C1 = 10
C2 = 2
position_list = [0, 1e-10, 2e-10]

marker = itertools.cycle(('o', 'v', 's', '8', 'p', '^', '<', '>'))
line = itertools.cycle(("-", "--", ":", "-."))

zeta_array = np.zeros([len(position_list), len(log10c_list)])
for j in range(0, len(log10c_list)):
    c = 10**log10c_list[j]
    c_H = 10**-pH
    c_Cl = c + c_H
    c_list = np.array([c_Cl, c, c_H])
    K_list = np.array([K_ads, 10**pKa])
    z_list = np.array([-1, 1, 1])
    v_list = np.array([False, True, True])
    sol = Solution_1plate(c_list, K_list, z_list, v_list, pH_effect=False, C_1=C1, C_2=C2)
    sol.solver_sigma()
    sol.solver_PB()
    interp = interp1d(sol.x, sol.psi, kind='cubic')
    for i in range(0, len(position_list)):
        zeta_array[i, j] = interp(position_list[i])

for i in range(0, len(position_list)):
    plt.plot(log10c_list, zeta_array[i, :]*1000, linestyle=next(line), label='$x_{\zeta}$=%.2f Ang' % (position_list[i]*1e10))
plt.xlabel('log$_{10}$(c) [M]')
plt.ylabel('Zeta potential [mV]')
plt.legend()
plt.text(-3, -40, 'p$K_{ads}$=%.2f\np$K_a$=%.2f' % (np.log10(K_ads), pKa))
plt.show()
