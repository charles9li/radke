import numpy as np
from scipy.constants import e, epsilon_0, k, N_A, R
from scipy.interpolate import interp1d
from scipy.integrate import quad, solve_bvp, trapz
from scipy.optimize import fsolve, root


class _Solution:
    """
    Parent class for Solution1Plate and Solution2Plate classes.

    To solve equations, a _Solution class is first instantiated with the
    desired parameters. The solve_equations method of the newly-instantiated
    _Solution object is called to solve relevant equations. A 3-pt quadratic
    continuation is employed to aid in convergence.

    Runs continuation on the following parameters:
        1. concentrations                           (c_list)
        2. dissociation constants                   (K_list)
        3. capacitance of inner Helmholtz region    (C1)
        4. capacitance of outer Helmholtz region    (C2)
        5. plate separation distance, for 2plate    (D)

    Parameters
    ----------
    c_list : (N,) array_like
             A 1-D array of concentrations. The anion concentration must be the
             first element.
    K_list : (N-1,) array_like
             A 1-D array of dissociation constants for each cation.
    z_list : (N,) array_like
             A 1-D array of valences of each ion.
    v_list : (N,) array_like
             A 1-D array of booleans that indicated whether or not each ion
             adsorbs the to surface.
    D : float
        Separation distance between two mica plates. Only used for the
        Solution2Plate class.
    C1 : float, optional
         Capacitance of the inner Helmholtz region. If a cation(s) is
         specified, then this value will be overridden.
    C2 : float, optional
         Capacitance of the outer Helmholtz region. If a cation(s) is
         specified, then this value will be overridden.
    pH_effect : bool, optional
                If true, then H+ complexes in the 0-plane.
    pH : float, optional
         pH of the bulk solution. Only used when pH_effect is set to True.
    pKa : float, optional
          pKa of the mica surface. Only used when pH_effect is set to True.
    T : float, optional
        Temperature of the system in Kelvin.
    L : float, optional
        Number of surface sites on mica basal plane
    eps_r_1 : float, optional
              Relative permitivitty of inner Helmholtz region.
    eps_r_2 : float, optional
              Relative permitivitty of outer Helmholtz region.
    eps_r_bulk : float, optional
                 Relative permitivitty of bulk solution.
    cation : str or (N,) array_like, optional
             Specifies cations in the system. If this parameters is not None,
             then dissociation constants and capacitances are computed.
              
    """

    # Dictionary of dissociation constants
    K_dict = {"Li": 10 ** 0.3,
              "Na": 2.5 * 10 ** 0.3,
              "K": 10 ** 2.8,
              "Cs": 10 ** 3}

    # Dictionary of crystallographic radii
    R_cr_dict = {"Li":  60e-12,
                 "Na":  98e-12,
                 "K":  133e-12,
                 "Cs": 169e-12}

    # Dictionary of hydrated radii
    R_hyd_dict = {"Li": 3.82e-10,
                  "Na": 3.52e-10,
                  "K":  3.31e-10,
                  "Cs": 3.29e-10}

    #########################
    # Initializing Instance #
    #########################

    def __init__(self, c_list, K_list, z_list, v_list, D, C1=0.5, C2=0.5,
                 pH_effect=True, pH=5.8, pKa=5.3, T=298, L=2e18, eps_r_1=6,
                 eps_r_2=30, eps_r_bulk=80, cation=None):
        """
        Initialize a _Solution instance.

        """
        self.pH_effect = pH_effect
        self.c_list = np.array(c_list)
        self.K_list = np.array(K_list)
        self._change_K_list(cation)
        self.z_list = np.array(z_list)
        self.v_list = np.array(v_list)
        self.D = D
        self.C1 = self._compute_C1(C1, cation, eps_r_1)
        self.C2 = self._compute_C2(C2, cation, eps_r_2)
        self.pH = pH
        self.pKa = pKa
        self.T = T
        self.L = L
        self.eps_r_1 = eps_r_1
        self.eps_r_2 = eps_r_2
        self.eps_r_bulk = eps_r_bulk
        self.eps_bulk = eps_r_bulk * epsilon_0
        self.cation = cation
        if pH_effect:
            self.c_list = np.append(self.c_list, 10 ** -pH)
            self.c_list[0] += 10 ** -pH
            self.K_list = np.append(self.K_list, 10 ** pKa)
            self.z_list = np.append(self.z_list, 1)
            self.v_list = np.append(self.v_list, True)

    def _compute_C1(self, C1, cation, eps_r_1):
        """
        Computes capacitance of inner Helmholtz region.

        If cation(s) are specified, then C1 is computed using crystallographic
        radius and the relative permitivitty of the inner Helmholtz region. If
        cation is None, then the specified value for C1 is used.

        """
        if cation is None:
            return C1
        eps_1 = eps_r_1*epsilon_0
        d1 = self._compute_d1(cation)
        return eps_1/d1

    def _compute_d1(self, cation):
        """
        Computes width of inner Helmholtz region. Uses the average
        crystallographic radius.

        """
        R_cr_avg = self._compute_R_cr_avg(cation)
        return R_cr_avg

    def _compute_C2(self, C2, cation, eps_r_2):
        """
        Computes capacitance of outer Helmholtz region.

        If cation(s) are specified, then C1 is computed using hydrated radius
        and the relative permitivitty of the outer Helmholtz region. If cation
        is None, then the specified value for C1 is used.

        """
        if cation is None:
            return C2
        eps_2 = eps_r_2*epsilon_0
        d2 = self._compute_d2(cation)
        return eps_2/d2

    def _compute_d2(self, cation):
        """
        Computes width of outer Helmholtz region. Uses the average hydrated
        radius.

        """
        R_hyd_avg = self._compute_R_hyd_avg(cation)
        return 2*R_hyd_avg

    def _compute_R_cr_avg(self, cation):
        """
        Computes the average crystallographic radius of the cation(s)
        specified.

        """
        if type(cation) is str:
            self._check_cation_R_cr(cation)
            return self.R_cr_dict[cation]
        else:
            R_cr_list = []
            for c in cation:
                self._check_cation_R_cr(c)
                R_cr_list += [self.R_cr_dict[c]]
            return np.mean(R_cr_list)

    def _compute_R_hyd_avg(self, cation):
        """
        Computes the average hydrated radius of the cation(s) specified.

        """
        if type(cation) is str:
            self._check_cation_R_hyd(cation)
            return self.R_hyd_dict[cation]
        else:
            R_hyd_list = []
            for c in cation:
                self._check_cation_R_hyd(c)
                R_hyd_list += [self.R_hyd_dict[c]]
            return np.mean(R_hyd_list)

    def _check_cation_R_cr(self, cation):
        """
        Raises AssertionError if value for crystallographic radius does not
        exist for specified cation.

        """
        message = "'" + cation + "' has no available value for " \
                                 "crystallographic radius value."
        assert cation in self.R_cr_dict.keys(), message

    def _check_cation_R_hyd(self, cation):
        """
        Raises AssertionError if value for hydrated radius does not exist for
        specified cation.

        """
        message = "'" + cation + "' has no available value for hydrated " \
                                 "radius."
        assert cation in self.R_hyd_dict.keys(), message

    ##########################
    # Dissociation Constants #
    ##########################

    def _change_K_list(self, cation):
        if cation is not None:
            if type(cation) is str:
                self.K_list[0] = self.K_dict[cation]
            else:
                for i in range(len(cation)):
                    self.K_list[i] = self.K_dict[cation[i]]

    def _check_cation_K_ads(self, cation):
        message = "'" + cation + "' has no available value for adsorption " \
                                 "equilibrium constant."
        assert cation in self.K_dict.keys(), message

    def compute_kappa(self, rho_list):
        return np.sqrt(e ** 2 * np.sum(self.z_list**2*rho_list) / (self.eps_bulk * k * self.T))

    def solve_equations(self, solution=None):
        """Solves equations to find potential profiles and surface charge
        densities.

        Uses a continuation method for concentration, adsorption constant,
        and capacitance parameters.

        """

        # Prints header containing parameter information
        self._print_header()

        # Initialize starting point for continuation
        self._create_c_list_init(solution)
        self._create_K_list_init(solution)
        self._create_C1_init(solution)
        self._create_C2_init(solution)
        self._create_D_init(solution)

        # Initialize guess
        guess = self._create_guess(solution)
        print("guess created")

        # Run continuation for each parameter
        guess = self._continuation('c_list', guess, 0.01, True)
        guess = self._continuation('K_list', guess, 0.01, True)
        guess = self._continuation('C1', guess, 0.1, False)
        guess = self._continuation('C2', guess, 0.1, False)
        guess = self._continuation('D', guess, 1e-9, False)
        print()

        self._solver(guess, self.c_list, self.K_list, self.C1, self.C2, self.D, get_values=True)
        self.guess = guess

    def _continuation(self, parameter_str, guess, step_size, log):
        guess_list = [guess]
        while not self._is_parameter_done(parameter_str):
            self._increment_parameter(parameter_str, step_size, log)
            guess_next = self._guess_next(guess_list)
            guess = self._solver_init(guess_next)
            guess_list = self._add_guess(guess_list, guess)
        print("{0} continuation done".format(parameter_str))
        return guess

    def _print_header(self):
        self._print_header_pH()
        self._print_header_T_L()
        self._print_header_eps()
        self._print_header_c()
        self._print_header_K()
        self._print_header_z()
        self._print_header_v()
        self._print_header_D_C()

    def _print_header_pH(self):
        string = "pH_effect : {0}".format(self.pH_effect)
        if self.pH_effect:
            string += ", pH = {0:.1f}, pKa = {1:.1f}".format(self.pH, self.pKa)
        print(string)

    def _print_header_T_L(self):
        print("T = {0:.0f} K, L = {1:.1e} 1/m2".format(self.T, self.L))

    def _print_header_eps(self):
        string = "eps_r_1 = {0:.0f}".format(self.eps_r_1)
        string += ", eps_r_2 = {0:.0f}".format(self.eps_r_2)
        string += ", eps_r_bulk = {0:.0f}".format(self.eps_r_bulk)
        print(string)

    def _print_header_c(self):
        string = "c_list : {0} M".format(self.c_list)
        print(string)

    def _print_header_K(self):
        string = "K_list : {0}".format(self.K_list)

    def _print_header_z(self):
        string = "z_list : {0}".format(self.z_list)
        print(string)

    def _print_header_v(self):
        string = "v_list : {0}".format(self.v_list)
        print(string)

    def _print_header_D_C(self):
        string = "D = {0:.1f} nm".format(self.D*1e9)
        string += ", C1 = {0:.2f} F/m2".format(self.C1)
        string += ", C2 = {0:.2f} F/m2".format(self.C2)
        print(string)
        print("-"*len(string))

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
            if self.pH_effect:
                self._c_list_init = np.array([0.1] * (len(self.c_list[self.v_list])-1))
                self._c_list_init = np.append(self._c_list_init, 10**-5.8)
            else:
                self._c_list_init = np.array([1e-3] * len(self.c_list[self.v_list]))
            self._c_list_init = np.append(np.sum(self._c_list_init), self._c_list_init)
        else:
            self._c_list_init = solution.c_list
        self._c_list_index = 1

    def _create_K_list_init(self, solution):
        if solution is None:
            if self.pH_effect:
                self._K_list_init = np.array([1e2] * (len(self.K_list)-1))
                self._K_list_init = np.append(self._K_list_init, 10**5.3)
            else:
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


class Solution1Plate(_Solution):

    def __init__(self, c_list, K_list, z_list, v_list, C1=0.5, C2=0.5,
                 pH=5.8, pKa=5.3, pH_effect=True, T=298, L=2e18, eps_r_bulk=80,
                 cation=None):
        super().__init__(c_list, K_list, z_list, v_list, 10e-9, C1=C1, C2=C2,
                         pH=pH, pKa=pKa, pH_effect=pH_effect, T=T, L=L, eps_r_bulk=eps_r_bulk,
                         cation=cation)

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):
        def equations(guess):

            psi_d, SM_list = guess[0], guess[1:]

            # Compute sigma_d
            c_bulk_sum = 1000*np.sum(c_list)
            c_d_sum = 1000*np.sum(c_list*np.exp(-self.z_list*e*psi_d/(k*self.T)))
            sigma_d = -psi_d/abs(psi_d)*np.sqrt(2 * R * self.T * self.eps_bulk * (c_d_sum - c_bulk_sum))

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
                num_cat = len(self.c_list) - 2
                psi_d_poly = np.poly1d([-0.06682361, 0.16257741,
                                        -0.14451976, 0.05710006,
                                        -0.00943752])
                SM_poly = np.poly1d([-7.19062004e+17, 1.83332403e+18,
                                     -1.78817799e+18, 9.06035496e+17,
                                     2.97346101e+17])
                SH_poly = np.poly1d([5.89869591e+17, -1.50684919e+18,
                                     1.47781243e+18, -7.65105121e+17,
                                     1.66578958e+18])
                guess = np.array([SM_poly(num_cat * 0.1)] * num_cat)
                guess = np.append(psi_d_poly(num_cat * 0.1), guess)
                guess = np.append(guess, SH_poly(num_cat * 0.1))
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


class Solution2Plate(_Solution):

    def _solver(self, guess, c_list, K_list, C1, C2, D, get_values=False):

        rho_list = 1000*N_A*c_list

        sigma_d_guess = guess[0]
        D_guess = guess[1]
        sol_guess = guess[2]

        size = 50
        x_guess = np.linspace(0, D, size)
        if sol_guess is None:
            kappa = self.compute_kappa(rho_list)
            psi_guess = -sigma_d_guess / (self.eps_bulk * kappa) * np.exp(-kappa * x_guess[0:size // 2])
            psi_guess = np.concatenate((psi_guess, list(reversed(psi_guess))))
            dpsi_guess = sigma_d_guess / self.eps_bulk * np.exp(-kappa * x_guess[0:size // 2])
            dpsi_guess = np.concatenate((dpsi_guess, list(reversed(dpsi_guess))))
        elif type(sol_guess) is tuple:
            psi_guess = np.zeros(len(x_guess))
            dpsi_guess = np.zeros(len(x_guess))
            psi_guess_list = [sol_guess[i](x_guess/D*D_guess[i])[0] for i in range(len(D_guess))]
            dpsi_guess_list = [sol_guess[i](x_guess/D*D_guess[i])[1] for i in range(len(D_guess))]
            for i in range(len(x_guess)):
                x = [j for j in range(len(psi_guess_list))]
                y = [psi_guess_list[j][i] for j in range(len(psi_guess_list))]
                dy = [dpsi_guess_list[j][i] for j in range(len(psi_guess_list))]
                psi_guess[i] = np.poly1d(np.polyfit(x, y, len(psi_guess_list) - 1))(-1)
                dpsi_guess[i] = np.poly1d(np.polyfit(x, dy, len(dpsi_guess_list) - 1))(-1)
        else:
            psi_guess = sol_guess(x_guess/D*D_guess)[0]
            dpsi_guess = sol_guess(x_guess/D*D_guess)[1]
        psi_guess = np.vstack((psi_guess, dpsi_guess))

        def solve_ode(sigma_d):

            def fun(x, psi):
                rho = np.zeros(len(x))
                for i in range(len(rho_list)):
                    rho_bulk = rho_list[i]
                    z = self.z_list[i]
                    rho += z*rho_bulk*np.exp(-z*e*psi[0]/(k*self.T))

                dpsi = psi[1]
                d2psi = -e / self.eps_bulk * rho
                return np.vstack((dpsi, d2psi))

            def bc(psia, psib):
                return np.array([psia[1] - sigma_d / self.eps_bulk, psib[1] + psia[1]])

            res = solve_bvp(fun, bc, x_guess, psi_guess)
            sol = res.sol

            if get_values:
                psi_beta = sol(0)[0] - sigma_d/C2
                sigma_0, sigma_beta = equations(sigma_d, psi_beta)
                self.sigma_0 = sigma_0
                self.sigma_beta = sigma_beta
                self.sigma_d = sigma_d
                self.sol = sol
                self.psi_m = sol(D/2)[0]
                self.P = self._compute_P()
            else:
                return sol

        def equations(sigma_d, psi_beta):
            S = sigma_d/e

            if self.pH_effect:
                self._assert_no_pH_effect()
            else:
                sigma_0 = -e*self.L
                SM_list = K_list*c_list[self.v_list]*S*np.exp(-self.z_list[self.v_list]*e*psi_beta/(k*self.T))
                SH = 0
                sigma_beta = sum(SM_list)*e
            return sigma_0, sigma_beta

        def objective(sigma_d):
            sol = solve_ode(sigma_d)
            psi_beta = sol(0)[0] - sigma_d/C2
            sigma_0, sigma_beta = equations(sigma_d, psi_beta)
            return sigma_0 + sigma_beta + sigma_d

        if get_values:
            solve_ode(sigma_d_guess)
        else:
            sigma_d = fsolve(objective, sigma_d_guess)[0]
            sol = solve_ode(sigma_d)
            return sigma_d, D, sol

    def _create_guess(self, solution):
        if solution is None:
            num_cat = len(self.c_list) - 1
            sigma_d_reg = np.poly1d([0.23274843, -0.07449596, 0.01438773])
            sigma_d = sigma_d_reg(num_cat * 1e-2)
            guess = (sigma_d, 10e-9, None)
            guess = self._solver_init(guess)
        else:
            guess = solution.guess
        return guess

    def _continuation_linear(self, guess_list):
        guess = [0, 0, 0]
        lin_reg = np.poly1d(np.polyfit([0, 1],
                                       [guess_list[0][0], guess_list[1][0]],
                                       1))
        guess[0] = lin_reg(-1)
        guess[1] = (guess_list[0][1], guess_list[1][1])
        guess[2] = (guess_list[0][2], guess_list[1][2])
        return guess

    def _continuation_quadratic(self, guess_list):
        guess = [0, 0, 0]
        quad_reg = np.poly1d(np.polyfit([0, 1, 2],
                                        [guess_list[0][0], guess_list[1][0], guess_list[2][0]],
                                        2))
        guess[0] = quad_reg(-1)
        guess[1] = (guess_list[0][1], guess_list[1][1], guess_list[2][1])
        guess[2] = (guess_list[0][2], guess_list[1][2], guess_list[2][2])
        return guess

    def _create_c_list_init(self, solution):
        if solution is None:
            self._c_list_init = np.array([1e-2] * len(self.c_list[self.v_list]))
            self._c_list_init = np.append(np.sum(self._c_list_init), self._c_list_init)
        else:
            self._c_list_init = solution.c_list
        self._c_list_index = 1

    def _compute_P(self):
        rho_list = self.c_list*N_A*1000
        rho_m = sum(rho_list*np.exp(-self.z_list*e*self.psi_m/(k*self.T)))
        rho_bulk = sum(rho_list)
        return k*self.T*(rho_m-rho_bulk)

    def psi_profile(self, x):
        return self.sol(x)[0]


class Force:

    R_cr_dict = {"Li": 60e-12,
                 "Na": 98e-12,
                 "K": 133e-12,
                 "Cs": 169e-12}
    R_hyd_dict = {"Li": 3.82e-10,
                  "Na": 3.52e-10,
                  "K": 3.31e-10,
                  "Cs": 3.29e-10}

    def __init__(self, c_list, K_list, z_list, v_list, C1=0.5, C2=0.5,
                 pH=5.8, pKa=5.3, pH_effect=True, T=298, L=2e18, eps_r_1=6,
                 eps_r_2=30, eps_r_bulk=80, A=2.2e-20, cation=None):

        def compute_sol(D, solution):
            sol = Solution2Plate(c_list, K_list, z_list, v_list, D, C1=C1, C2=C2,
                                 pH=pH, pKa=pKa, pH_effect=pH_effect, T=T, L=L,
                                 eps_r_1=eps_r_1, eps_r_2=eps_r_2, eps_r_bulk=eps_r_bulk,
                                 cation=cation)
            sol.solve_equations(solution)
            return sol

        self._compute_sol = compute_sol
        self.R_cr = self.R_cr_dict[cation]
        self.R_hyd = self.R_hyd_dict[cation]
        self.A = A

    def compute_F(self, D_start):
        self.compute_W_osmotic(D_start)
        W_vdw_list = np.zeros(len(self.D_list))
        for i in range(len(W_vdw_list)):
            W_vdw_list[i] = -self.A/(12*np.pi*(self.D_list[i]+2*self.R_cr+4*self.R_hyd)**2)
        self.W_vdw_list = W_vdw_list
        self.W_list = self.W_osmotic_list + W_vdw_list
        self.F_list = 2*np.pi*self.W_list

    def compute_W_osmotic(self, D_start):
        D_list = np.array([])
        P_list = np.array([])
        sol_prev = None
        W_prev = -1
        W_curr = self._integrate_P(D_list, P_list)
        D = D_start

        while abs(W_curr-W_prev) > 1e-10 or len(D_list) < 3:
            print(D*1e9)
            sol = self._compute_sol(D, sol_prev)
            sol_prev = sol
            D_list = np.append(D_list, D)
            P_list = np.append(P_list, sol.P)
            W_prev = W_curr
            W_curr = self._integrate_P(D_list, P_list)
            D = self._increment_D(D)

        W_osmotic_list = np.zeros(len(D_list) - 1)
        P_interp = interp1d(D_list, P_list, kind='cubic')
        for i in range(len(W_osmotic_list)):
            W_osmotic = quad(P_interp, D_list[i], D_list[-1])[0]
            W_osmotic_list[i] = W_osmotic

        self.D_list = D_list[:-1]
        self.W_osmotic_list = W_osmotic_list

    @staticmethod
    def _integrate_P(D_list, P_list):
        if len(D_list) <= 1:
            return 1
        elif len(D_list) <= 3:
            return trapz(P_list, x=D_list)
        else:
            P_interp = interp1d(D_list, P_list, kind='cubic')
            W = quad(P_interp, D_list[0], D_list[-1])[0]
            return W

    @staticmethod
    def _increment_D(D):
        if D < 5e-9:
            return D + 0.1e-9
        else:
            return D + 1e-9
