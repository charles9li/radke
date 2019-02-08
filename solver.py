import numpy as np
import warnings
from constants import *
from scipy.constants import *
from scipy.optimize import root, minimize
from scipy.integrate import odeint, simps, solve_bvp


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

	# Solve Poisson-Boltzmann to get potential and ion distributions
	def solver_PB(self):

		# Checks to see if solver_sigma method has been successfully called
		if not self.solver_sigma_complete:
			self.solver_sigma()

		c_list = self.c_list
		z_list = self.z_list

		# Convert bulk concentrations to number density
		rho_list = 1000*N_A*c_list

		# Calculated values
		kappa = np.sqrt(e**2*np.sum(rho_list)/(epsilon_0*eps_bulk*k*T))

		def solver(sigma_d, x_end, psi_guess=None):

			def fun(x, psi):
				d2psi_dx2 = np.zeros(len(psi[0]))
				for i in range(len(rho_list)):
					d2psi_dx2 += z_list[i]*rho_list[i]*np.exp(-z_list[i]*e*psi[0]/(k*T))
				d2psi_dx2 *= -e/(epsilon_0*eps_bulk)
				dpsi_dx = psi[1]
				return np.vstack((dpsi_dx, d2psi_dx2))

			def bc(psia, psib):
				return np.array([psia[1]-sigma_d/(eps_0*eps_bulk), psib[0]])

			size = 50
			x_dist = np.linspace(0, x_end, size)

			if psi_guess is None:
				psi_guess = -sigma_d/(eps_0*eps_bulk*kappa)*np.exp(-kappa*x_dist)
				dpsi_guess = sigma_d/(eps_0*eps_bulk)*np.exp(-kappa*x_dist)
				psi_guess = np.vstack((psi_guess, dpsi_guess))

			res = solve_bvp(fun, bc, x_dist, psi_guess)
			psi = res.sol(x_dist)[0]
			dpsi = res.sol(x_dist)[1]
			return np.vstack((x_dist, psi, dpsi))

		# Sets overflow errors to warning
		np.seterr(over='warn')
		warnings.filterwarnings('error')

		# Increases x span until furthest psi value is close enough to 0	
		x_end = 1e-9
		sol = solver(self.sigma_d, x_end)
		self.x = sol[0]
		self.psi = sol[1]
		psi_guess = sol[1:]
		while abs((self.psi[0]-self.psi_d)/self.psi_d) > 1e-4:
			x_end += 1e-9
			sol = solver(self.sigma_d, x_end, psi_guess)
			self.x = sol[0]
			self.psi = sol[1]
			psi_guess = sol[1:]

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

		# Solve charge density and site balance equations to find sigma and beta values
	def solver_sigma_PB(self):

		c_list = self.c_list
		K_list = self.K_list
		z_list = self.z_list
		v_list = self.v_list

		# Solve equations until psi_m converges
		def solver(psi_m):

			def equations(X0, c_list, K_list, get_values=False):

				psi_d, SM_list = X0[0], X0[1:]

				# Compute sigma_d
				c_d_sum = 1000*np.sum(c_list*np.exp(-z_list*e*psi_d/(k*T)))
				c_m_sum = 1000*np.sum(c_list*np.exp(-z_list*e*psi_m/(k*T)))
				sigma_d = -psi_d/abs(psi_d)*np.sqrt(2*R*T*eps_bulk*epsilon_0*(c_d_sum - c_m_sum))

				# Compute psi_beta and number of free sites S
				psi_beta = psi_d - sigma_d/self.C_2
				S = L - np.sum(SM_list)

				# Create adsoprtion equations
				if self.pH_effect:
					sigma_0 = -e*(L-SM_list[-1])
					psi_0 = psi_beta + sigma_0/self.C_1
					SM_objective = (K_list[:-1] - SM_list[:-1]/(S*(c_list[:-1]*np.exp(-z_list[:-1]*e*psi_beta/(k*T)))[v_list[:-1]]))/K_list[:-1]
					SH_objective = (K_list[-1] - SM_list[-1]/(S*(c_list[-1]*np.exp(-z_list[-1]*e*psi_0/(k*T)))))/K_list[-1]
					SM_objective = np.append(SM_objective, SH_objective)
				else:
					sigma_0 = -e*L
					psi_0 = psi_beta + sigma_0/C_1
					SM_objective = (K_list - SM_list/(S*(c_list*np.exp(-z_list*e*psi_beta/(k*T)))[v_all]))/K_list

				# Create total charge equation
				if self.pH_effect:
					sigma_objective = (sigma_0 + sigma_d + e*np.sum((z_list[v_list]*SM_list)[:-1]))/sigma_0
				else:
					sigma_objective = (sigma_0 + sigma_d + e*np.sum(z_list[v_list*SM_list]))/sigma_0

				# Create soution attributes or return equations
				if get_values:
					self.SM_list = SM_list
					self.psi_0 = psi_0
					self.psi_beta = psi_beta
					self.psi_d = psi_d
					self.sigma_0 = sigma_0
					self.sigma_beta_list = e*SM_list
					self.sigma_d = sigma_d
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
				c_list_prev = np.array([0.1, 0.05, 0.05])
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

			# Calculated values
			rho_list = c_list*1000*N_A # bulk number density of Na+	(m^-3)
			kappa = np.sqrt(np.sum((e*z_list)**2*rho_list)/(epsilon_0*eps_bulk*k*T))

			# # Poisson-Boltzmann equation
			# def fun(psi, x):
			# 	c_m_sum = np.sum(c_list*1000*np.exp(-z_list*e*psi_m/(k*T)))
			# 	c_sum = np.sum(c_list*1000*np.exp(-z_list*e*psi/(k*T)))
			# 	dpsidx = np.sqrt(2*R*T/(epsilon_0*eps_bulk)*(c_sum-c_m_sum))
			# 	return dpsidx

			# size = 100
			# x_dist = np.linspace(0, D/2, size)
			# sol = odeint(fun, self.psi_d, x_dist)

			# Poisson-Boltzmann equation
			def fun(psi, x):
				psi1, psi2 = psi
				rho_sum = np.sum(z_list*rho_list*np.exp(-z_list*e*psi1/(k*T)))
				dpsi1dx = psi2
				dspi2dx = -e/(epsilon_0*eps_bulk)*rho_sum
				return [dpsi1dx, dspi2dx]

			size = 100
			x_dist = np.linspace(0, self.D/2, size)
			psi0 = [self.psi_d, self.sigma_d/(epsilon_0*eps_bulk)]
			sol = odeint(fun, psi0, x_dist)

			self.x = x_dist
			self.psi = sol[:,0]
			self.psi_m = psi_m

		# Finds initial condition for solver
		def objective(psi_m):
			solver(psi_m)
			print(psi_m)
			print(self.psi[-1])
			return (psi_m - self.psi[-1])**2

		# root(objective, 0)
		minimize(objective, 0, bounds=((-np.inf, 0),))

		# Indicates that the solver_sigma_PB method has been successfully run
		self.solver_sigma_PB_complete = True
















# pH = 5.8
# pKa = 5.3
# c_H = 10**-pH
# c_1 = 0.15000000000000002*0.67976
# c_2 = (1-0.15000000000000002)*0.67976
# c_Cl = c_H + c_1 + c_2
# c_list = [c_Cl, c_1, c_2]
# K_list = [10**2.4, 100]
# z_list = [-1, 1, 1]
# v_list = [False, True, True]
# D = 10e-9

# pH = 5.8
# pKa = 5.3
# c_H = 10**-pH
# c = 1e-3
# c_Cl = c + c_H
# c_list = [c_Cl, c, c_H]
# K_list = [10**2.4, 10**5.3]
# z_list = [-1, 1, 1]
# v_list = [False, True, True]
# D = 10e-9

# sol = Solution_2plate(c_list, K_list, z_list, v_list, D, pKa=pKa)
# sol.solver_sigma_PB()
# print(sol.psi_0)
# print(sol.psi_beta)
# print(sol.psi_d)
# print(sol.sigma_0)
# print(sol.sigma_beta_list)
# print(sol.sigma_d)
