import matplotlib.pyplot as plt
from scipy.constants import e, k
from data_pashley import *
from class_Pashley import *

# Figure 1
# psi_d_list_fig1 = np.ones(len(pashley_fig1_data))
# for data in pashley_fig1_data:
#     c = data.c
data = pashley_fig1_data[0]
c = data.c
D = data.D
FR = data.FR
plt.scatter(D, FR*1e-3)
sol = PashleySolution(c, -0.05, np.linspace(D[0], D[-1])*1e-9)
sol.compute_FR()
plt.semilogy(sol.D_list*1e9, sol.FR_list*1e3)
plt.show()
