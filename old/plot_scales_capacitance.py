import matplotlib.pyplot as plt
from old.class_Scales import *

# p*KM = 3.0, vary C1
from old.line_marker import *
plt.figure('pKM3_varyC1')
C1_list = np.array([0.5, 1, 5, 10])
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
fontsize = 14
plt.rc('font', size=fontsize)
plt.rc('axes', labelsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
plt.text(-2.5, -120, 'pH = %.1f\np*$K_M$ = %.1f\n$C_2$ = %i $\mu$C/cm$^2$' % (sol.pH, sol.pKM, sol.C2*100))
plt.legend(title='$C_1$\n[$\mu$C/cm$^2$]')

# p*KM = 3.0, vary C2
from old.line_marker import *
plt.figure('pKM3_varyC2')
C2_list = np.array([0.5, 1, 5, 10])
for C2 in C2_list:
    c_list = np.logspace(-5, -1)
    zeta_list = np.ones(len(c_list))
    i = 0
    for c in c_list:
        sol = Solution_Scales(c, 3.0, C2=C2)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(lines), label='%i' % (C2*100))
fontsize = 14
plt.rc('font', size=fontsize)
plt.rc('axes', labelsize=fontsize)
plt.rc('xtick', labelsize=fontsize)
plt.rc('ytick', labelsize=fontsize)
plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
plt.text(-2.5, -120, 'pH = %.1f\np*$K_M$ = %.1f\n$C_2$ = %i $\mu$C/cm$^2$' % (sol.pH, sol.pKM, sol.C1*100))
plt.legend(title='$C_2$\n[$\mu$C/cm$^2$]')

# Show all plots
plt.show()
