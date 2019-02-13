import warnings
import matplotlib.pyplot as plt
from scipy.constants import e, k, R, epsilon_0
from scipy.optimize import root
from data_scales import *


class Solution_Scales:

    solver_sigma_complete = False

    def __init__(self, c, pKM, pH=5.8, pKa=5.6, C1=10, C2=2, Ns=2e18, T=293, eps=79):
        self.c = c
        self.pH = pH
        self.cH = 10**-pH
        self.cCl = self.c + self.cH
        self.KM = 10**-pKM
        self.Ka = 10**-pKa
        self.C1 = C1
        self.C2 = C2
        self.Ns = Ns
        self.T = T
        self.eps = eps

    def solver_sigma(self):

        def equations(X0, c, KM, pH, get_values=False):

            # Unpack variables
            psi_d, SM, SH = X0

            # Calculate cH and cCl
            cH = 10**-pH
            cCl = c + cH

            # Calculate sigma_d
            c_bulk_sum = 1000*np.sum([cCl, c, cH])
            c_d_sum = 1000*((c + cH)*np.exp(-e*psi_d/(k*self.T)) + cCl*np.exp(e*psi_d/(k*self.T)))
            sigma_d = -psi_d/abs(psi_d)*np.sqrt(2*R*self.T*self.eps*epsilon_0*(c_d_sum - c_bulk_sum))

            # Compute quantities
            psi_beta = psi_d - sigma_d/self.C2
            S = self.Ns - SM - SH
            sigma_beta = e*SM
            sigma_0 = -e*(S + SM)
            psi_0 = psi_beta + sigma_0/self.C1

            # # Equations
            # S = self.Ns - SM - SH
            # sigma_0 = -e*(SM + S)
            # sigma_beta = e*SM
            # sigma_d = e*S
            # psi_beta = psi_d - sigma_d/self.C2
            # psi_0 = psi_beta + sigma_0/self.C1

            # Return objective or create solution attributes
            if get_values:
                self.SM = SM
                self.SH = SH
                self.S = S
                self.sigma_0 = sigma_0
                self.sigma_beta = sigma_beta
                self.sigma_d = sigma_d
                self.psi_0 = psi_0
                self.psi_beta = psi_beta
                self.psi_d = psi_d
            else:
                sigma_objective = (sigma_0 + sigma_beta + sigma_d)/sigma_0
                KM_objective = (KM - SM*cH/(SH*c)*np.exp(e*(psi_beta-psi_0)/(k*self.T)))/KM
                Ka_objective = (self.Ka - S*cH/SH*np.exp(-e*psi_0/(k*self.T)))/self.Ka
                return np.array([sigma_objective, KM_objective, Ka_objective])

        # Set overflow to trigger warning
        np.seterr(over='warn')
        warnings.filterwarnings('error')

        # Helper function to create new guesses
        def guess_create(equations, guess, c, KM, pH):
            root_func = lambda X0: equations(X0, c, KM, pH)
            solution = root(root_func, guess, method='lm', tol=1e-10)
            return solution.x

        # Helper function to take log mean of K
        def log_mean(K1, K2):
            return 10**np.mean([np.log10(K1), np.log10(K2)])

        # Initialize guess and starting c and KM values
        X0 = np.array([-0.005973597351999846, 1.8921500128751857e+18, 8.052532531933648e+16])
        c_prev = 0.1
        KM_prev = 10**-2.9
        pH_prev = 5.8
        X0 = guess_create(equations, X0, c_prev, KM_prev, pH_prev)
        c_curr = np.mean([self.c, c_prev])
        KM_curr = log_mean(self.KM, KM_prev)
        pH_curr = np.mean([self.pH, pH_prev])

        # Iterate through c, KM, and pH values to update guess until convergence
        success = False
        while not success:
            try:
                X0 = guess_create(equations, X0, self.c, self.KM, self.pH)
                equations(X0, self.c, self.KM, self.pH, get_values=True)
                success = True
            except Warning:
                try:
                    X0 = guess_create(equations, X0, c_curr, KM_curr, pH_curr)
                    c_prev = c_curr
                    KM_prev = KM_curr
                    pH_prev = pH_curr
                    c_curr = np.mean([self.c, c_curr])
                    KM_curr = log_mean(self.KM, KM_curr)
                    pH_curr = np.mean([self.pH, pH_curr])
                except Warning:
                    c_curr = np.mean([c_prev, c_curr])
                    KM_curr = log_mean(KM_prev, KM_curr)
                    pH_curr = np.mean([pH_prev, pH_curr])

        # Indications that the solver sigma method has been successfully run
        self.solver_sigma_complete = True


# Set to true for plot
fig1 = False
fig2 = True
fig10 = True
fig11 = True

# Figure 2
if fig2:
    plt.figure('plot_scales_fig2')
    plt.scatter(pH_data_fig2, zeta_data_fig2)
    fontsize = 14
    pH_list = np.linspace(2, 12)
    zeta_list = np.ones(len(pH_list))
    i = 0
    for pH in pH_list:
        sol = Solution_Scales(1e-3, 3.2, pH=pH)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    plt.plot(pH_list, zeta_list*1000)
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.title('Scales, Figure 2')
    plt.xlabel('pH', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)

# Figure 10
if fig10:
    from line_marker import *
    plt.figure('plot_scales_fig10')
    salts = ['CsCl', 'KCl']
    pKM = dict(CsCl=2.8, KCl=3.2)
    for salt in salts:
        c_data = c_data_fig1[salt]
        zeta_data = zeta_data_fig1[salt]
        plt.scatter(np.log10(c_data), zeta_data, marker=next(markers), label=salt)
        c_list = np.logspace(-4.5, -2)
        zeta_list = np.ones(len(c_list))
        i = 0
        for c in c_list:
            sol = Solution_Scales(c, pKM[salt])
            sol.solver_sigma()
            zeta_list[i] = sol.psi_d
            i += 1
        plt.plot(np.log10(c_list), zeta_list*1000, linestyle=next(lines))
    fontsize = 14
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.title('Scales, Figure 10')
    plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.text(-4.5, -55, 'pH=5.8')
    plt.legend()

# Figure 11
if fig11:
    from line_marker import *
    plt.figure('plot_scales_fig11')
    salts = ['NaCl', 'LiCl']
    for salt in salts:
        c_data = c_data_fig1[salt]
        zeta_data = zeta_data_fig1[salt]
        plt.scatter(np.log10(c_data), zeta_data, marker=next(markers), label=salt)
    c_list = np.logspace(-4, -2, 10)
    zeta_list = np.ones(len(c_list))
    i = 0
    for c in c_list:
        sol = Solution_Scales(c, 5.5)
        sol.solver_sigma()
        zeta_list[i] = sol.psi_d
        i += 1
    fontsize = 14
    plt.rc('font', size=fontsize)
    plt.rc('axes', labelsize=fontsize)
    plt.rc('xtick', labelsize=fontsize)
    plt.rc('ytick', labelsize=fontsize)
    plt.plot(np.log10(c_list), zeta_list*1000, color='k', linestyle=':')
    plt.title('Scales, Figure 11')
    plt.xlabel('log$_{10}$(c) [mol/L]', fontsize=fontsize)
    plt.ylabel('Zeta potential [mV]', fontsize=fontsize)
    plt.text(-4, -78, 'pH=5.8')
    plt.legend()

# Show all figures
plt.show()
