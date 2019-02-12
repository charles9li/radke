import numpy as np
import warnings
from scipy.constants import e, k
from scipy.optimize import root

class Solution_Scales:

    solver_sigma_complete = False

    def __init__(self, c, pKM, pH=5.8, pKa=5.6, C1=10, C2=2, Ns=2e18, T=293):
        self.c = c
        self.cH = 10**-pH
        self.cCl = self.c + self.cH
        self.KM = 10**-pKM
        self.Ka = 10**-pKa
        self.C1 = C1
        self.C2 = C2
        self.Ns = Ns
        self.T = T

    def solver_sigma(self):

        def equations(X0, c, KM, get_values=False):

            psi_d, SM, SH = X0

            # Equations
            S = self.Ns - SM - SH
            sigma_0 = -e*(SM + S)
            sigma_beta = e*SM
            sigma_d = e*S
            psi_beta = psi_d - sigma_d/self.C2
            psi_0 = psi_beta + sigma_0/self.C1

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
                KM_objective = self.KM - SM*self.cH/(SH*self.c)*np.exp(e*(psi_beta-psi_0)/(k*self.T))
                Ka_objective = self.Ka - S*self.cH/SH*np.exp(-e*psi_0/(k*self.T))
                return np.array([sigma_objective, KM_objective, Ka_objective])

        # Set overflow to trigger warning
        np.seterr(over='warn')
        warnings.filterwarnings('error')

        # Helper function to create new guesses
        def guess_create(equations, guess, c, KM):
            root_func = lambda X0: equations(X0, c, KM)
            solution = root(root_func, guess, method='lm', tol=1e-10)
            return solution.x

        # Helper function to take log mean of K
        def log_mean(K1, K2):
            return 10**np.mean([np.log10(K1), np.log10(K2)])

        # Initialize guess and starting c and KM values
        X0 = np.array([-0.01, 0.2*self.Ns, 0.2*self.Ns])
        c_prev = 0.1
        KM_prev = 10**-3
        X0 = guess_create(equations, X0, c_prev, KM_prev)
        print(X0)
        c_curr = np.mean([self.c, c_prev])
        KM_curr = log_mean(self.KM, KM_prev)

        # Iterate through K and c values to update guess until convergence
        success = False
        while not success:
            try:
                X0 = guess_create(equations, X0, self.c, self.KM)
                equations(X0, self.c, self.KM, get_values=True)
                success = True
            except Warning:
                print('Warning')
                try:
                    X0 = guess_create(equations, X0, c_curr, KM_curr)
                    c_prev = c_curr
                    KM_prev = KM_curr
                    c_curr = np.mean([self.c, c_curr])
                    KM_curr = log_mean(self.KM, KM_curr)
                except Warning:
                    c_curr = np.mean([c_prev, c_curr])
                    KM_curr = log_mean([KM_prev, KM_curr])

        # Indications that the solver sigma method has been successfully run
        self.solver_sigma_complete = True


sol = Solution_Scales(.1, 2.9)
sol.solver_sigma()
print(sol.psi_0)
print(sol.psi_beta)
print(sol.psi_d)
print(sol.sigma_0)
print(sol.sigma_beta)
print(sol.sigma_d)
