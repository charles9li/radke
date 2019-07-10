import matplotlib.pyplot as plt
import itertools
from solver_old import *

# Parameters
log10c_list = np.linspace(-5, 0, 10)
K_ads = 10**3.151010098187411
K_ads = 100
pH = 5.8
pKa = 5.3
C1_list = [0.2, 0.5, 1, 2, 5, 10]
C2_list = [0.5, 1, 2]

marker = itertools.cycle(('o', 'v', 's', '8', 'p', '^', '<', '>'))
line = itertools.cycle(("-", "--", ":", "-."))

C2 = 0.2
plt.figure(1)
for C1 in C1_list:
    zeta_list = np.zeros(len(log10c_list))
    i = 0
    for log10c in log10c_list:
        c = 10**log10c
        c_H = 10**-pH
        c_Cl = c + c_H
        c_list = np.array([c_Cl, c, c_H])
        K_list = np.array([K_ads, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list, K_list, z_list, v_list, pH_effect=False, C_1=C1, C_2=C2)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(
        log10c_list, zeta_list*1000, linestyle=next(line), label='$C_1$=%d, $C_2$=%d' % (C1*100, C2*100))
plt.xlabel('log$_{10}$(c) [M]')
plt.ylabel('Zeta potential [mV]')
plt.legend(title='Capacitances [$\mu$F/cm$^2$]')
plt.text(-3, -40, 'p$K_{ads}$=%.2f\np$K_a$=%.2f' % (np.log10(K_ads), pKa))

C1 = .2
plt.figure(2)
for C2 in C2_list:
    zeta_list = np.zeros(len(log10c_list))
    i = 0
    for log10c in log10c_list:
        c = 10**log10c
        c_H = 10**-pH
        c_Cl = c + c_H
        c_list = np.array([c_Cl, c, c_H])
        K_list = np.array([K_ads, 10**pKa])
        z_list = np.array([-1, 1, 1])
        v_list = np.array([False, True, True])
        sol = Solution_1plate(c_list, K_list, z_list, v_list, pH_effect=False, C_1=C1, C_2=C2)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(
        log10c_list, zeta_list*1000, linestyle=next(line), label='$C_1$=%d, $C_2$=%d' % (C1*100, C2*100))
plt.xlabel('log$_{10}$(c) [M]')
plt.ylabel('Zeta potential [mV]')
plt.legend(title='Capacitances [$\mu$F/cm$^2$]')
plt.text(-3, -40, 'p$K_{ads}$=%.2f\np$K_a$=%.2f' % (np.log10(K_ads), pKa))

plt.show()
