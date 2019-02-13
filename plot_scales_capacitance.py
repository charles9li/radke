import matplotlib.pyplot as plt
from class_Scales import *

# p*KM = 3.0, vary C1
from line_marker import *
plt.figure('pKM3_varyC1')
C1_list = np.linspace(2, 10, 5)
C1_list = np.array([0.1, 0.25, 0.5, 1, 5, 10])
for C1 in C1_list:
    c_list = np.logspace(-5, -1)
    zeta_list = np.ones(len(c_list))
    i = 0
    for c in c_list:
        sol = Solution_Scales(c, 3.0, C1=C1)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(lines), label='%i' % (C1*100))
fontsize = 15
plt.rc('font', size=fontsize)
plt.rc('axes', labelsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
plt.text(-3.8, -35, 'pH = 5.8\np*$K_M$ = 3.0\n$C_2$ = %i $\mu$C/cm$^2$' % (sol.C2*100))
plt.legend(title='$C_1$\n[$\mu$C/cm$^2$]')

plt.show()
