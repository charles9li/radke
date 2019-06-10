import numpy as np
import warnings
from scipy.constants import e, k, epsilon_0, R
from scipy.optimize import root, minimize
from scipy.integrate import odeint, simps, solve_bvp
np.seterr(over='warn')
warnings.filterwarnings('error')


class Solution:

    L = 2e18

    def __init__(self, c_list, K_list, z_list, v_list, D, pH=5.8, pKa=5.3, pH_effect=True, C1=0.5, C2=0.5, T=298):
        self.c_list = np.array(c_list)
        self.K_list = np.array(K_list)
        self.z_list = np.array(z_list)
        self.v_list = np.array(v_list)
        self.D = D
        self.pH = pH
        self.pKa = pKa
        self.pH_effect = pH_effect
        self.C1 = C1
        self.C2 = C2
        self.T = T
        if pH_effect:
            self.c_list = np.append(self.c_list, 10**-pH)
            self.K_list = np.append(self.K_list, 10**pKa)
            self.z_list = np.append(self.z_list, 1)
            self.v_list = np.append(self.v_list, True)

    def solve_equations(self, c_list_init=None, K_list_init=None,
                        C1_init=None, C2_init=None, D_init=None, guess=None):
        """Solves equations to find potential profiles and surface charge
        densities.

        Uses a continuation method for concentration, adsorption constant,
        and capacitance parameters.

        """

        # Initialize starting point for continuation
        self._create_c_list_init(c_list_init)
        self._create_K_list_init(K_list_init)
        self._create_C1_init(C1_init)
        self._create_C2_init(C2_init)
        self._create_D_init(D_init)

        # Initialize guess
        guess = self._create_guess(guess)

        # Run continuation for each parameter
        guess = self._continuation('c_list', self._log_mean, guess)
        guess = self._continuation('K_list', self._log_mean, guess)
        guess = self._continuation('C1', self._mean, guess)
        guess = self._continuation('C2', self._mean, guess)
        guess = self._continuation('D', self._mean, guess)

        self._solver(guess, self.c_list, self.K_list, self.C1, self.C2, self.D, get_values=True)

    def _continuation(self, parameter_str, average, guess):
        while True:
            try:
                guess_full = self._solver(guess, *self._parameter_full(parameter_str))
                guess_half_1 = self._solver(guess, *self._parameter_half(parameter_str, average))
                guess_half_2 = self._solver(guess, *self._parameter_full(parameter_str))

                if self._is_guess_converged(guess_full, guess_half_2):
                    return guess_full
                else:
                    self._parameter_change_curr(parameter_str, average)
                    warnings.warn("{} continuation. Guess did not converge".format(parameter_str), Warning)

            except Warning:
                try:
                    guess_full = self._solver(guess, *self._parameter_full(parameter_str))
                    guess_half_1 = self._solver(guess, *self._parameter_half(parameter_str, average))
                    guess_half_2 = self._solver(guess, *self._parameter_full(parameter_str))

                    if not self._is_guess_converged(guess_full, guess_half_2):
                        warnings.warn("{} continuation. Guess did not converge".format(parameter_str), Warning)

                    self._parameter_change_init(parameter_str)
                    self._parameter_target(parameter_str)
                except Warning:
                    self._parameter_change_curr(parameter_str, average)

    def _parameter_full(self, parameter_str):
        c_list = self._c_list_init
        K_list = self._K_list_init
        C1 = self._C1_init
        C2 = self._C2_init
        D = self._D_init
        if parameter_str == 'c_list':
            c_list = self._c_list_curr
        elif parameter_str == 'K_list':
            K_list = self._K_list_curr
        elif parameter_str == 'C1':
            C1 = self._C1_curr
        elif parameter_str == 'C2':
            C2 = self._C2_curr
        elif parameter_str == 'D':
            D = self._D_curr
        return c_list, K_list, C1, C2, D

    def _parameter_half(self, parameter_str, average):
        c_list = self._c_list_init
        K_list = self._K_list_init
        C1 = self._C1_init
        C2 = self._C2_init
        D = self._D_init
        if parameter_str == 'c_list':
            c_list = average(c_list, self._c_list_curr)
            c_list[0] = np.sum(c_list[1:])
        elif parameter_str == 'K_list':
            K_list = average(K_list, self._K_list_curr)
        elif parameter_str == 'C1':
            C1 = average(C1, self._C1_curr)
        elif parameter_str == 'C2':
            C2 = average(C2, self._C2_curr)
        elif parameter_str == 'D':
            D = average(D, self._D_curr)
        return c_list, K_list, C1, C2, D

    def _parameter_target(self, parameter_str):
        if parameter_str == 'c_list':
            self._c_list_curr = self.c_list
        elif parameter_str == 'K_list':
            self._K_list_curr = self.K_list
        elif parameter_str == 'C1':
            self._C1_curr = self.C1
        elif parameter_str == 'C2':
            self._C2_curr = self.C2
        elif parameter_str == 'D':
            self._D_curr = self.D

    def _parameter_change_curr(self, parameter_str, average):
        if parameter_str == 'c_list':
            self._c_list_curr = average(self._c_list_init, self._c_list_curr)
            self._c_list_curr[0] = np.sum(self._c_list_curr[1:])
        elif parameter_str == 'K_list':
            self._K_list_curr = average(self._K_list_init, self._K_list_curr)
        elif parameter_str == 'C1':
            self._C1_curr = average(self._C1_init, self._C1_curr)
        elif parameter_str == 'C2':
            self._C2_curr = average(self._C2_init, self._C2_curr)
        elif parameter_str == 'D':
            self._D_curr = average(self._D_init, self._D_curr)

    def _parameter_change_init(self, parameter_str):
        if parameter_str == 'c_list':
            self._c_list_init = self._c_list_curr
        elif parameter_str == 'K_list':
            self._K_list_init = self._K_list_curr
        elif parameter_str == 'C1':
            self._C1_init = self._C1_curr
        elif parameter_str == 'C2':
            self._C2_init = self._C2_curr
        elif parameter_str == 'D':
            self._D_init = self._D_curr

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):
        pass

    #
    # Initializes starting c and K values for continuation
    #

    def _create_c_list_init(self, c_list_init):
        if c_list_init is None:
            self._c_list_init = np.array([1e-3] * len(self.c_list[self.v_list]))
            self._c_list_init = np.append([np.sum(self._c_list_init)], self._c_list_init)
        else:
            self._c_list_init = c_list_init
        self._c_list_curr = self.c_list[self.v_list]

    def _create_K_list_init(self, K_list_init):
        if K_list_init is None:
            self._K_list_init = np.array([1e2] * len(self.K_list))
        else:
            self._K_list_init = K_list_init
        self._K_list_curr = self.K_list

    def _create_C1_init(self, C1_init):
        if C1_init is None:
            self._C1_init = 0.5
        else:
            self._C1_init = C1_init
        self._C1_curr = self.C1

    def _create_C2_init(self, C2_init):
        if C2_init is None:
            self._C2_init = 0.5
        else:
            self._C2_init = C2_init
        self._C2_curr = self.C2

    def _create_D_init(self, D_init):
        if D_init is None:
            self._D_init = 10e-9
        else:
            self._D_init = D_init
        self._D_curr = self.D

    def _create_guess(self, guess):
        return guess

    #
    # Utility methods
    #

    @staticmethod
    def _log_mean(x, y):
        return 10**Solution._mean(np.log10(x), np.log10(y))

    @staticmethod
    def _mean(x, y):
        return (x + y) / 2

    @staticmethod
    def _is_guess_converged(guess_1, guess_2):
        return False


class Solution1Plate(Solution):

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):
        def equations(guess):

            psi_d, SM_list = guess[0], guess[1:]

            # Compute sigma_d
            c_bulk_sum = 1000*np.sum(c_list)
            c_d_sum = 1000*np.sum(c_list*np.exp(-self.z_list*e*psi_d/(k*self.T)))
            sigma_d = -psi_d/abs(psi_d)*np.sqrt(2*R*self.T*79*epsilon_0*(c_d_sum - c_bulk_sum))

            # Compute psi_beta and number of free sites S
            psi_beta = psi_d - sigma_d/C2
            S = self.L - np.sum(SM_list)

            # Create adsorption equations
            if self.pH_effect:
                sigma_0 = -e*(self.L-SM_list[-1])
                psi_0 = psi_beta + sigma_0/C1
                SM_objective = (K_list[:-1] - SM_list[:-1]/(S*(c_list[:-1]*np.exp(-self.z_list[:-1]*e*psi_beta/(k*self.T)))[self.v_list[:-1]]))/K_list[:-1]
                SH_objective = (K_list[-1] - SM_list[-1]/(S*(c_list[-1]*np.exp(-self.z_list[-1]*e*psi_0/(k*self.T)))))/K_list[-1]
                SM_objective = np.append(SM_objective, SH_objective)
            else:
                sigma_0 = -e*self.L
                psi_0 = psi_beta + sigma_0/C1
                SM_objective = (K_list - SM_list/(S*(c_list*np.exp(-self.z_list*e*psi_beta/(k*self.T)))[self.v_list]))/K_list

            # Create total charge equation
            if self.pH_effect:
                sigma_objective = (sigma_0 + sigma_d + e*np.sum((self.z_list[self.v_list]*SM_list)[:-1]))/sigma_0
            else:
                sigma_objective = (sigma_0 + sigma_d + e*np.sum(self.z_list[self.v_list]*SM_list))/sigma_0

            # Create solution attributes or return equations
            if get_values:
                self.SM_list = SM_list
                self.psi_0 = psi_0
                self.psi_beta = psi_beta
                self.psi_d = psi_d
                self.sigma_0 = sigma_0
                self.sigma_beta_list = e*SM_list
                self.sigma_d = sigma_d
                self.frac_list = self.SM_list/2e18
                if self.pH_effect:
                    self.SM_list = SM_list[:-1]
                    self.SH = SM_list[-1]
                    self.sigma_beta_list = e*SM_list[:-1]
            else:
                return np.append(SM_objective, sigma_objective)

        if get_values:
            equations(guess)
        else:
            solution = root(equations, guess, method='lm', tol=1e-10)
            return solution.x

    def _create_guess(self, guess):
        if guess is None:
            return np.array([-0.10727440871184632, 1.9076140644838958e18])
        return guess

    @staticmethod
    def _is_guess_converged(guess_1, guess_2):
        i = 0
        for i in range(len(guess_1)):
            item_1 = guess_1[i]
            item_2 = guess_2[i]
            if abs(item_1 - item_2) > 1e-3:
                return False
        return True


class Solution2Plate(Solution):
    pass


c_list = [1e-3, 1e-3]
K_list = [100]
z_list = [-1, 1]
v_list = [False, True]
sol = Solution1Plate(c_list, K_list, z_list, v_list, 10e-9, pH_effect=False)
sol.solve_equations()
