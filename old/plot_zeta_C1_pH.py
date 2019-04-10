import matplotlib.pyplot as plt
import itertools
from solvers import *

# Parameters
c_list = 10**np.linspace(-5, -1, 10)
pH = 5.8
K_ads = 100
pKa = 5.6
C1_list = np.array([.1, .2, .5, 1, 5, 10])
C2 = .2

marker = itertools.cycle(('o', 'v', 's', '8', 'p', '^', '<', '>'))
line = itertools.cycle(("-", "--", ":", "-."))

i = 0
for C1 in C1_list:
    zeta_list = np.zeros(len(c_list))
    j = 0
    for c in c_list:
        sol = solver_1plate_1cation_ads_pH(c, pH, K_ads, pKa, C_1=C1, C_2=C2)
        zeta_list[j] = sol.psi[2]
        j += 1
    plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(line), label='$C_1$ = %i $\mu$F/cm$^2$' % (100*C1))
    i += 1

plt.xlabel('log$_{10}$(c) [M]')
plt.ylabel('Zeta potential [mV]')
plt.legend()
plt.show()
