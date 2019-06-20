import numpy as np
from scipy.constants import e, k, epsilon_0, R, N_A
from scipy.optimize import root, fsolve
from scipy.integrate import odeint, simps, solve_bvp


class Solution:

    SM_list = None
    guess = None

    def __init__(self, c_list, K_list, z_list, v_list, D, C1=0.5, C2=0.5,
                 pH=5.8, pKa=5.3, pH_effect=True, T=298, L=2e18, eps_r=80):
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
        self.L = L
        self.eps = eps_r * epsilon_0
        if pH_effect:
            self.c_list = np.append(self.c_list, 10**-pH)
            self.K_list = np.append(self.K_list, 10**pKa)
            self.z_list = np.append(self.z_list, 1)
            self.v_list = np.append(self.v_list, True)

    def compute_kappa(self, c_list):
        return np.sqrt(e**2*np.sum(self.z_list**2*self.c_list)/(self.eps*k*self.T))

    def solve_equations(self, solution=None):
        """Solves equations to find potential profiles and surface charge
        densities.

        Uses a continuation method for concentration, adsorption constant,
        and capacitance parameters.

        """

        # Initialize starting point for continuation
        self._create_c_list_init(solution)
        self._create_K_list_init(solution)
        self._create_C1_init(solution)
        self._create_C2_init(solution)
        self._create_D_init(solution)

        # Initialize guess
        guess = self._create_guess(solution)

        # Run continuation for each parameter
        guess = self._continuation('c_list', guess, 0.01, True)
        guess = self._continuation('K_list', guess, 0.01, True)
        guess = self._continuation('C1', guess, 0.1, False)
        guess = self._continuation('C2', guess, 0.1, False)
        guess = self._continuation('D', guess, 1e-9, False)

        self._solver(guess, self.c_list, self.K_list, self.C1, self.C2, self.D, get_values=True)
        self.guess = guess

    def _continuation(self, parameter_str, guess, step_size, log):
        guess_list = [guess]
        while not self._is_parameter_done(parameter_str):
            self._increment_parameter(parameter_str, step_size, log)
            guess_next = self._guess_next(guess_list)
            guess = self._solver_init(guess_next)
            guess_list = self._add_guess(guess_list, guess)
        return guess

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):
        pass

    def _solver_init(self, guess):
        return self._solver(guess, self._c_list_init, self._K_list_init,
                            self._C1_init, self._C2_init, self._D_init)

    #
    # Used for continuation
    #

    @staticmethod
    def _add_guess(guess_list, guess):
        if len(guess_list) < 3:
            return [guess, *guess_list]
        else:
            guess_list[1:] = guess_list[0:2]
            guess_list[0] = guess
            return guess_list

    def _guess_next(self, guess_list):
        if len(guess_list) == 1:
            return self._solver_init(guess_list[0])
        elif len(guess_list) == 2:
            return self._continuation_linear(guess_list)
        return self._continuation_quadratic(guess_list)

    def _continuation_linear(self, guess_list):
        pass

    def _continuation_quadratic(self, guess_list):
        pass

    def _is_parameter_done(self, parameter_str):
        return_bool = False
        if parameter_str == 'c_list':
            return_bool = self._c_list_index >= len(self._c_list_init)
        elif parameter_str == 'K_list':
            return_bool = self._K_list_index >= len(self._K_list_init)
        elif parameter_str == 'C1':
            return_bool = self._C1_init == self.C1
        elif parameter_str == 'C2':
            return_bool = self._C2_init == self.C2
        elif parameter_str == 'D':
            return_bool = self._D_init == self.D
        if return_bool:
            print(parameter_str + " success")
        return return_bool

    #
    # Increments parameters
    #

    def _increment_parameter(self, parameter_str, step_size, log=False):
        if parameter_str == 'c_list':
            self._increment_c(step_size, log)
        elif parameter_str == 'K_list':
            self._increment_K(step_size, log)
        elif parameter_str == 'C1':
            self._C1_init = self._increment(self._C1_init, self.C1, step_size, log)
        elif parameter_str == 'C2':
            self._C2_init = self._increment(self._C2_init, self.C2, step_size, log)
        elif parameter_str == 'D':
            self._D_init = self._increment(self._D_init, self.D, step_size, log)

    def _increment_c(self, step_size, log):
        init = self._c_list_init[self._c_list_index]
        final = self.c_list[self._c_list_index]

        self._c_list_init[self._c_list_index] = self._increment(init, final, step_size, log)
        self._c_list_init[0] = sum(self._c_list_init[1:])

        init = self._c_list_init[self._c_list_index]
        if init == final:
            self._c_list_index += 1

    def _increment_K(self, step_size, log):
        init = self._K_list_init[self._K_list_index]
        final = self.K_list[self._K_list_index]

        self._K_list_init[self._K_list_index] = self._increment(init, final, step_size, log)

        init = self._K_list_init[self._K_list_index]
        if init == final:
            self._K_list_index += 1

    @staticmethod
    def _increment(init, final, step_size, log):
        if init == final:
            return init

        sign = (final - init) / abs(final - init)

        if log:
            init = 10**(np.log10(init) + sign * step_size)
        else:
            init = init + sign * step_size

        if sign * init > sign * final:
            return final
        else:
            return init

    #
    # Initializes starting c and K values for continuation
    #

    def _create_c_list_init(self, solution):
        if solution is None:
            self._c_list_init = np.array([1e-3] * len(self.c_list[self.v_list]))
            self._c_list_init = np.append(np.sum(self._c_list_init), self._c_list_init)
        else:
            self._c_list_init = solution.c_list
        self._c_list_index = 1

    def _create_K_list_init(self, solution):
        if solution is None:
            self._K_list_init = np.array([1e2] * len(self.K_list))
        else:
            self._K_list_init = solution.K_list
        self._K_list_index = 0

    def _create_C1_init(self, solution):
        if solution is None:
            self._C1_init = 0.5
        else:
            self._C1_init = solution.C1

    def _create_C2_init(self, solution):
        if solution is None:
            self._C2_init = 0.5
        else:
            self._C2_init = solution.C2

    def _create_D_init(self, solution):
        if solution is None:
            self._D_init = 10e-9
        else:
            self._D_init = solution.D

    def _create_guess(self, solution):
        return solution

    #
    # Utility methods
    #

    @staticmethod
    def _assert_no_pH_effect():
        assert False, "pH effect not implemented"


class Solution1Plate(Solution):

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):
        def equations(guess):

            psi_d, SM_list = guess[0], guess[1:]

            # Compute sigma_d
            c_bulk_sum = 1000*np.sum(c_list)
            c_d_sum = 1000*np.sum(c_list*np.exp(-self.z_list*e*psi_d/(k*self.T)))
            sigma_d = -psi_d/abs(psi_d)*np.sqrt(2*R*self.T*self.eps*(c_d_sum - c_bulk_sum))

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
                self.frac_list = self.SM_list/self.L
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

    def _continuation_linear(self, guess_list):
        guess = np.zeros(len(guess_list[0]))
        for i in range(len(guess)):
            lin_reg = np.poly1d(np.polyfit([0, 1],
                                           [guess_list[0][i], guess_list[1][i]],
                                           1))
            guess[i] = lin_reg(-1)
        return guess

    def _continuation_quadratic(self, guess_list):
        guess = np.zeros(len(guess_list[0]))
        for i in range(len(guess)):
            quad_reg = np.poly1d(np.polyfit([0, 1, 2],
                                 [guess_list[0][i], guess_list[1][i], guess_list[2][i]],
                                 2))
            guess[i] = quad_reg(-1)
        return guess

    def _create_guess(self, solution):
        if solution is None:
            if self.pH_effect:
                self._assert_no_pH_effect()
            else:
                num_cat = len(self.c_list) - 1
                psi_d_poly = np.poly1d([3.25590388e+09, -1.08057231e+08,
                                        1.41687447e+06, -9.49887189e+03,
                                        3.72769964e+01, -1.36321197e-01])
                SM_poly = np.poly1d([-1.52703963e+34,  6.65368881e+32,
                                     -1.21054552e+31,  1.19304587e+29,
                                     -6.90913117e+26,  2.38311837e+24,
                                     -4.72002436e+21,  4.82749852e+18])
            guess = np.array([SM_poly(num_cat * 1e-3)] * num_cat)
            guess = np.append(psi_d_poly(num_cat * 1e-3), guess)
            guess = self._solver(guess, self._c_list_init, self._K_list_init,
                                 self._C1_init, self._C2_init, self._D_init)
        else:
            guess = solution.guess
        return guess


class Solution2Plate(Solution):

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):

        rho_list = 1000*N_A*c_list

        sigma_d_guess = guess[0]
        D_guess = guess[1]
        sol_guess = guess[2]

        size = 50
        if sol_guess is None:
            x_guess = np.linspace(0, D, size)
            kappa = self.compute_kappa(c_list)
            psi_guess = -sigma_d_guess/(self.eps*kappa)*np.exp(-kappa*x_guess[0:size//2])
            psi_guess = np.concatenate((psi_guess, list(reversed(psi_guess))))
            dpsi_guess = sigma_d_guess/self.eps*np.exp(-kappa*x_guess[0:size//2])
            dpsi_guess = np.concatenate((dpsi_guess, list(reversed(dpsi_guess))))
            psi_guess = np.vstack((psi_guess, dpsi_guess))
        else:
            sigma_d = guess[0]
            x_guess = np.linspace(0, D, size)
            psi_guess = guess[2]

        def solve_ode(guess):

            def fun(x, psi):
                rho = np.zeros(len(x))
                for i in range(len(rho_list)):
                    rho_bulk = rho_list[i]
                    z = self.z_list[i]
                    rho += rho_bulk*np.exp(-z*e*psi[0]/(k*self.T))

                dpsi = psi[1]
                d2psi = -e/self.eps*rho
                return np.vstack((dpsi, d2psi))

            def bc(psia, psib):
                return np.array([psia[1] - sigma_d/self.eps, psib[1] + psia[1]])

            res = solve_bvp(fun, bc, x_guess, psi_guess)
            return res.sol

        def equations(sigma_d, psi_beta):
            S = sigma_d/e

            if self.pH_effect:
                self._assert_no_pH_effect()
            else:
                sigma_0 = -e*self.L
                SM_list = K_list*c_list[self.v_list]*S*np.exp(-self.z_list[self.v_list]*e*psi_beta/(k*self.T))
                SH = 0
                sigma_beta = sum(SM_list)*e

            if get_values:
                self.sigma_0 = sigma_0
                self.SM_list = SM_list
                self.SH = SH
            else:
                return sigma_0, sigma_beta

        def objective(sigma_d):
            pass

        if get_values:
            pass
        else:
            sigma_d = fsolve(objective, sigma_d_guess)
            return sigma_d, D

    def _create_guess(self, solution):
        if solution is None:
            guess = (0, 10e-9, None)
        else:
            guess = solution.guess
        return guess

    def _guess_next(self, guess_list):
        pass
