import numpy as np
import warnings
from scipy.constants import e, k, epsilon_0
from scipy.optimize import root, minimize
from scipy.integrate import odeint, simps, solve_bvp
np.seterr(over='warn')
warnings.filterwarnings('error')


class Solution:

    def __init__(self, c_list, K_list, z_list, v_list, D, pH=5.8, pKa=5.3, pH_effect=True, C1=10, C2=2):
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

    def _continuation(self, parameter_str, average, guess):
        while True:
            try:
                guess_full = self._solver(guess, *self._parameter_full(parameter_str))
                guess_half_1 = self._solver(guess, *self._parameter_half(parameter_str))
                guess_half_2 = self._solver(guess, *self._parameter_full(parameter_str))
            except Warning:
                pass

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
        elif parameter_str == 'K_list':
            K_list = average(K_list, self._K_list_curr)
        elif parameter_str == 'C1':
            C1 = average(C1, self._C1_curr)
        elif parameter_str == 'C2':
            C2 = average(C2, self._C2_curr)
        elif parameter_str == 'D':
            D = average(D, self._D_curr)
        return c_list, K_list, C1, C2, D

    def _solver(self, guess, c_list, K_list, C1, C2, D):
        pass

    #
    # Initializes starting c and K values for continuation
    #

    def _create_c_list_init(self, c_list_init):
        if c_list_init is None:
            self._c_list_init = np.array([1e-3] * len(self.c_list[self.v_list]))
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
    # Continuation utility
    #

    def _create_argument(self):
        pass

    #
    # Utility methods
    #

    @staticmethod
    def _log_mean(x, y):
        return 10**Solution._mean(np.log10(x), np.log10(y))

    @staticmethod
    def _mean(x, y):
        return (x + y) / 2


class Solution1Plate(Solution):
    pass


class Solution2Plate(Solution):
    pass
