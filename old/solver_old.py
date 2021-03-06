import numpy as np
import warnings
from old.constants import *
from scipy.constants import e, k, epsilon_0
from scipy.optimize import root
from scipy.integrate import simps


class Solution_1plate:

    solver_sigma_complete = False
    solver_PB_complete = False
    bound_diffuse_complete = False

    def __init__(self, c_list, K_list, z_list, v_list, pH=5.8, pKa=5.3, pH_effect=True, C_1=10, C_2=2):
        self.c_list = c_list
        self.K_list = K_list
        self.z_list = z_list
        self.v_list = v_list
        self.pH = pH
        self.pKa = pKa
        self.pH_effect = pH_effect
        self.C_1 = C_1
        self.C_2 = C_2
        self.psi_0 = None
        self.psi_beta = None
        self.psi_d = None
        self.SM_list = None
        self.SH = None
        if pH_effect:
            self.c_list = np.append(self.c_list, 10**-pH)
            self.K_list = np.append(self.K_list, 10**pKa)
            self.z_list = np.append(self.z_list, 1)
            self.v_list = np.append(self.v_list, True)

    # Solve charge density and site balance equations to find sigma and beta values
    def solver_sigma(self):

        c_list = self.c_list
        K_list = self.K_list
        z_list = self.z_list
        v_list = self.v_list

        def equations(X0, c_list, K_list, get_values=False):

            psi_d, SM_list = X0[0], X0[1:]

            # Compute sigma_d
            c_bulk_sum = 1000*np.sum(c_list)
            c_d_sum = 1000*np.sum(c_list*np.exp(-z_list*e*psi_d/(k*T)))
            sigma_d = -psi_d/abs(psi_d)*np.sqrt(2*R*T*eps_bulk*epsilon_0*(c_d_sum - c_bulk_sum))

            # Compute psi_beta and number of free sites S
            psi_beta = psi_d - sigma_d/self.C_2
            S = L - np.sum(SM_list)

            # Create adsorption equations
            if self.pH_effect:
                sigma_0 = -e*(L-SM_list[-1])
                psi_0 = psi_beta + sigma_0/self.C_1
                SM_objective = (K_list[:-1] - SM_list[:-1]/(S*(c_list[:-1]*np.exp(-z_list[:-1]*e*psi_beta/(k*T)))[v_list[:-1]]))/K_list[:-1]
                SH_objective = (K_list[-1] - SM_list[-1]/(S*(c_list[-1]*np.exp(-z_list[-1]*e*psi_0/(k*T)))))/K_list[-1]
                SM_objective = np.append(SM_objective, SH_objective)
            else:
                sigma_0 = -e*L
                psi_0 = psi_beta + sigma_0/C_1
                SM_objective = (K_list - SM_list/(S*(c_list*np.exp(-z_list*e*psi_beta/(k*T)))[v_list]))/K_list

            # Create total charge equation
            if self.pH_effect:
                sigma_objective = (sigma_0 + sigma_d + e*np.sum((z_list[v_list]*SM_list)[:-1]))/sigma_0
            else:
                sigma_objective = (sigma_0 + sigma_d + e*np.sum(z_list[v_list]*SM_list))/sigma_0

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

        # Set overflow to trigger warning
        np.seterr(over='warn')
        warnings.filterwarnings('error')

        # Helper function to create new guesses
        def guess_create(equations, guess, c_list, K_list):
            root_func = lambda X0: equations(X0, c_list, K_list)
            solution = root(root_func, guess, method='lm', tol=1e-10)
            return solution.x

        # Helper function to take log mean of K lists
        def log_mean(K_list_1, K_list_2):
            return 10**np.mean([np.log10(K_list_1), np.log10(K_list_2)], axis=0)

        # Initialize guess and starting c and K values
        guess = np.append(-0.01, np.zeros(len(c_list[v_list])) + 0.2*L)
        if self.pH_effect:
            K_list_prev = np.append(np.ones(len(K_list[:-1])), 10**5.3)
            c_list_prev = np.append(0.1, 0.1/(len(c_list[:-2]))*np.ones(len(c_list[:-2])))
            c_list_prev = np.append(c_list_prev, 10**-5.8)
        else:
            K_list_prev = np.ones(len(c_list[v_list]))
            c_list_prev = np.append(0.1, 0.1/(len(c_list[:-1]))*np.ones(len(c_list[:-1])))
        guess = guess_create(equations, guess, c_list_prev, K_list_prev)
        c_list_curr = np.mean([c_list_prev, c_list], axis=0)
        K_list_curr = log_mean(K_list, K_list_prev)

        # Iterate through K and c values to update guess until convergence
        success = False
        while not success:
            try:
                X0 = guess_create(equations, guess, c_list, K_list)
                equations(X0, c_list, K_list, get_values=True)
                success = True
            except Warning:
                try:
                    guess = guess_create(equations, guess, c_list_curr, K_list_curr)
                    c_list_prev = c_list_curr
                    K_list_prev = K_list_curr
                    c_list_curr = np.mean([c_list, c_list_curr], axis=0)
                    K_list_curr = log_mean(K_list, K_list_curr)
                except Warning:
                    c_list_curr = np.mean([c_list_prev, c_list_curr], axis=0)
                    K_list_curr = log_mean(K_list_prev, K_list_curr)

        # Indicates that the solver_sigma method has been successfully run
        self.solver_sigma_complete = True

    # # Solve Poisson-Boltzmann to get potential and ion distributions
    # def solver_PB(self):
    #
    #     # Checks to see if solver_sigma method has been successfully called
    #     if not self.solver_sigma_complete:
    #         self.solver_sigma()
    #
    #     c_list = self.c_list
    #     z_list = self.z_list
    #
    #     # Convert bulk concentrations to number density
    #     rho_list = 1000*N_A*c_list
    #
    #     # Calculated values
    #     kappa = np.sqrt(e**2*np.sum(rho_list)/(epsilon_0*eps_bulk*k*T))
    #
    #     def solver(sigma_d, x_end, psi_guess=None):
    #
    #         def fun(x, psi):
    #             d2psi_dx2 = np.zeros(len(psi[0]))
    #             for i in range(len(rho_list)):
    #                 d2psi_dx2 += z_list[i]*rho_list[i]*np.exp(-z_list[i]*e*psi[0]/(k*T))
    #             d2psi_dx2 *= -e/(epsilon_0*eps_bulk)
    #             dpsi_dx = psi[1]
    #             return np.vstack((dpsi_dx, d2psi_dx2))
    #
    #         def bc(psia, psib):
    #             return np.array([psia[1]-sigma_d/(eps_0*eps_bulk), psib[0]])
    #
    #         size = 50
    #         x_dist = np.linspace(0, x_end, size)
    #
    #         if psi_guess is None:
    #             psi_guess = -sigma_d/(eps_0*eps_bulk*kappa)*np.exp(-kappa*x_dist)
    #             dpsi_guess = sigma_d/(eps_0*eps_bulk)*np.exp(-kappa*x_dist)
    #             psi_guess = np.vstack((psi_guess, dpsi_guess))
    #
    #         res = solve_bvp(fun, bc, x_dist, psi_guess)
    #         psi = res.sol(x_dist)[0]
    #         dpsi = res.sol(x_dist)[1]
    #         return np.vstack((x_dist, psi, dpsi))
    #
    #     # Sets overflow errors to warning
    #     np.seterr(over='warn')
    #     warnings.filterwarnings('error')
    #
    #     # Increases x span until furthest psi value is close enough to 0
    #     x_end = 1e-9
    #     sol = solver(self.sigma_d, x_end)
    #     self.x = sol[0]
    #     self.psi = sol[1]
    #     psi_guess = sol[1:]
    #     while abs((self.psi[0]-self.psi_d)/self.psi_d) > 1e-4:
    #         x_end += 1e-9
    #         sol = solver(self.sigma_d, x_end, psi_guess)
    #         self.x = sol[0]
    #         self.psi = sol[1]
    #         psi_guess = sol[1:]
    #
    #     # Creates an array of ion number density profiles
    #     self.rho_list = [rho_list[i]*np.exp(-z_list[i]*e*self.psi/(k*T)) for i in range(len(rho_list))]
    #
    #     # Converts each element in all of the number density profiles to float type
    #     self.rho_list = [np.array([float(r) for r in rho]) for rho in self.rho_list]
    #
    #     # Indicates that the solver_PB method has been successfully run
    #     self.solver_PB_complete = True

    # Solve Poisson-Boltzmann to get potential and ion distributions
    def solver_PB(self):

        # Checks to see if solver_sigma method has been successfully called
        if not self.solver_sigma_complete:
            self.solver_sigma()

        c_list = self.c_list
        z_list = self.z_list

        # Convert bulk concentrations to number density
        rho_list = 1000*N_A*c_list

        # Differential equation
        def fun(x, psi):
            return np.sqrt(2*k*T/(eps_bulk*epsilon_0)*np.sum(rho_list*(np.exp(-z_list*e*psi/(k*T))-1)))

        # Compute solution at next step
        def y_next(f, t, y, h):
            k1 = h*f(t,     y)
            k2 = h*f(t+h/2, y+k1/2)
            k3 = h*f(t+h/2, y+k2/2)
            k4 = h*f(t+h,   y+k3)
            return y + 1/6*(k1+2*k2+2*k3+k4)

        # Initiate x and psi solutions
        x = [0]
        psi = [self.psi_d]

        # Default step size
        h_default = 1e-9

        # Solve ODE
        while psi[-1] < -1e-6:
            h = h_default
            y_next_2h = y_next(fun, x[-1], psi[-1], 2*h)
            y_next_h1 = y_next(fun, x[-1], psi[-1], h)
            y_next_h2 = y_next(fun, x[-1]+h, y_next_h1, h)
            while abs(y_next_2h - y_next_h2) > 1e-8:
                h = h/2
                y_next_2h = y_next_h1
                y_next_h1 = y_next(fun, x[-1], psi[-1], h)
                y_next_h2 = y_next(fun, x[-1]+h, y_next_h1, h)
            x += [x[-1]+h]
            psi += [y_next_h1]

        # Store x and psi as instance variables
        self.x = np.array(x)
        self.psi = np.array(psi)

        # Creates an array of ion number density profiles
        self.rho_list = [rho_list[i]*np.exp(-z_list[i]*e*self.psi/(k*T)) for i in range(len(rho_list))]

        # Converts each element in all of the number density profiles to float type
        self.rho_list = [np.array([float(r) for r in rho]) for rho in self.rho_list]

        # Indicates that the solver_PB method has been successfully run
        self.solver_PB_complete = True


    # Compute surface density of each ion bound in the diffuse layer
    def bound_diffuse(self):

        # Checks to see if solver_PB was called successfully
        if not self.solver_PB_complete:
            self.solver_PB()

        # Converts bulk concentrations from mol/L to m^-3
        rho_bulk_list = 1000*N_A*self.c_list

        # Computes surface density of each ion bound in the diffuse layer
        self.bound_diffuse_list = [simps(self.rho_list[i]-rho_bulk_list[i], self.x) for i in range(len(self.rho_list))]

        # Indicates that the bound_diffuse method has been successfully run
        self.bound_diffuse_complete = True






class Solution_2plate:

    solver_sigma_PB_complete = False

    def __init__(self, c_list, K_list, z_list, v_list, D, pH=5.8, pKa=5.3, pH_effect=True, C_1=10, C_2=10):
        self.c_list = c_list
        self.K_list = K_list
        self.z_list = z_list
        self.v_list = v_list
        self.D = D
        self.pH = pH
        self.pKa = pKa
        self.pH_effect = pH_effect
        self.C_1 = C_1
        self.C_2 = C_2
        if pH_effect:
            self.c_list = np.append(self.c_list, 10**-pH)
            self.K_list = np.append(self.K_list, 10**pKa)
            self.z_list = np.append(self.z_list, 1)
            self.v_list = np.append(self.v_list, True)

    def solve(self):
        pass



c_list = np.array([1e-3, 1e-3])
K_list = np.array([100])
z_list = np.array([-1, 1])
v_list = np.array([False, True])
sol = Solution_1plate(c_list, K_list, z_list, v_list, pH_effect=False, C_1=0.5, C_2=0.5)
sol.solver_sigma()
