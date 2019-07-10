"""
SOLVERS
=======

This module contains a number of solvers for cases including 1- and 2-plate and
1- and 2-cation.

"""

# numpy and scipy

# for detecting overflow errors during iteration

# solution classes, constants, and helper function
from old.helper import *
from old.constants import *
from scipy.integrate import solve_bvp
from old.class_Solution import *
from scipy.interpolate import interp1d
import warnings
from scipy.optimize import fsolve

def solver_1plate_1cation_ads_pH(c, pH, K_ads, pKa, **kwargs):
	
	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	x_end = 10e-9

	def solver(sigma_d, c, x_end, psi_guess=None):

		# Calculated values
		rho_Na_bulk = c*1000*N_A 	# bulk number density of Na+	(m^-3)
		rho_H_bulk = 10**(-pH)*1000*N_A
		rho_Cl_bulk = rho_Na_bulk+rho_H_bulk 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(2*e**2*(rho_Na_bulk+rho_H_bulk+rho_Cl_bulk)/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_H = rho_H_bulk*np.exp(-e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_Na+rho_H+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
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
		rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_H = rho_H_bulk*np.exp(-e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_1plate(x_dist, -(x_2+x_1), -x_2)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		sigma_0, sigma_beta, psi_0 = sigma_0_beta_calc_1cation_pH(c, pH, K_ads, pKa, L, C_1, C_2, psi_beta)
		psi = append_1plate(psi, psi_0, psi_beta)
		dpsi = append_1plate(dpsi, 0, 0)
		rho_Na = append_1plate(rho_Na, 0, 0)
		rho_H = append_1plate(rho_H, 0, 0)
		rho_Cl = append_1plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_Na, rho_H, rho_Cl, dpsi))
		# sigma_0, sigma_beta = sigma_0_beta_calc_1cation_pH(c, pH, K_ads, pKa, L, psi_beta)
		return Solution_1cation_pH(solution, sigma_0, sigma_beta, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c, x_end, psi_guess=None):
		solution = solver(sigma_d, c, x_end, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	# Sets overflow errors to warning
	np.seterr(over='warn')
	warnings.filterwarnings('error')

	# Print progress header
	# print_1plate_1cation_ads_pH_header(c, pH, K_ads, pKa, C_1, C_2)

	# Initialize guess
	guess_exist = False
	while not guess_exist:
		try:
			sigma_d_guess = fsolve(lambda sigma_d: objective(sigma_d, c, x_end), 0)
			solution = solver(sigma_d_guess, c, x_end)
			psi_guess = solution.psi_guess_1plate()
			guess_exist = True
			# print_1plate_ads(x_end, sigma_d_guess)
		except Warning:
			x_end = x_end/2

	# Iterates until dpsi/dx = 0 as x -> inf
	while not solution_end(solution):
		x_end_prev = x_end
		x_end += 10e-9
		good = False
		while not good:
			try:
				sigma_d_guess = fsolve(
					lambda sigma_d: objective(sigma_d, c, x_end, psi_guess),
					sigma_d_guess)
				solution = solver(sigma_d_guess, c, x_end, psi_guess)
				# print_1plate_ads(x_end, sigma_d_guess)
				psi_guess = solution.psi_guess_1plate()
				good = True
			except Warning:
				x_end = np.mean([x_end_prev, x_end])
	return solution





def solver_1plate_2cation_ads(c_1, c_2, K_ads_1, K_ads_2, **kwargs):

	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter_string_list:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	x_end = 10e-9

	# Total cation concentration
	c_tot = c_1+c_2

	# Computes cation 1 and cation 2 fractions in bulk
	frac_1 = c_1/c_tot
	frac_2 = c_2/c_tot

	# Solves for potential and ion distrubutions in diffuse region
	def solver(sigma_d, c_1, c_2, x_end, psi_guess=None):

		# Calculated values
		rho_1_bulk = c_1*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_2_bulk = c_2*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_Cl_bulk = c_tot*1000*N_A 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(e**2*(rho_1_bulk+rho_2_bulk+rho_Cl_bulk)/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_1+z_Li*rho_2+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
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
		rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_1plate(x_dist, -(x_2+x_1), -x_2)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		psi_0 = psi_0_calc(psi_beta, sigma_0, C_1)
		psi = append_1plate(psi, psi_0, psi_beta)
		dpsi = append_1plate(dpsi, 0, 0)
		rho_1 = append_1plate(rho_1, 0, 0)
		rho_2 = append_1plate(rho_2, 0, 0)
		rho_Cl = append_1plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_1, rho_2, rho_Cl, dpsi))
		sigma_beta, sigma_beta_1, sigma_beta_2 = sigma_beta_calc_2cation(c_1, c_2, K_ads_1, K_ads_2, L, psi_beta)
		return Solution_2cation(solution, sigma_0, sigma_beta, sigma_beta_1, sigma_beta_2, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c_1, c_2, x_end, psi_guess=None):
		solution = solver(sigma_d, c_1, c_2, x_end, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	# Sets overflow errors to warning
	np.seterr(over='warn')
	warnings.filterwarnings('error')

	# Print progress header
	print_1plate_2cation_ads_header(c_1, c_2, K_ads_1, K_ads_2, C_1, C_2)

	# Initialize guess
	guess_exist = False
	while not guess_exist:
		try:
			sigma_d_guess = fsolve(lambda sigma_d: objective(sigma_d, c_1, c_2, x_end), 0)
			solution = solver(sigma_d_guess, c_1, c_2, x_end)
			psi_guess = solution.psi_guess_1plate()
			guess_exist = True
			print_1plate_ads(x_end, sigma_d_guess)
		except Warning:
			x_end = x_end/2

	# Iterates until dpsi/dx = 0 as x -> inf
	while not solution_end(solution):
		x_end_prev = x_end
		x_end += 10e-9
		good = False
		while not good:
			try:
				sigma_d_guess = fsolve(
					lambda sigma_d: objective(sigma_d, c_1, c_2, x_end, psi_guess),
					sigma_d_guess)
				solution = solver(sigma_d_guess, c_1, c_2, x_end, psi_guess)
				print_1plate_ads(x_end, sigma_d_guess)
				psi_guess = solution.psi_guess_1plate()
				good = True
			except Warning:
				x_end = np.mean([x_end_prev, x_end])
	return solution





def solver_1plate_2cation_ads_pH(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, eng, **kwargs):

	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter_string_list:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	interval = 1e-9

	x_end = interval

	# H+ concentration
	c_H = 10**(-pH)

	# Total cation concentration
	c_tot = c_1+c_2+c_H

	# Computes cation 1 and cation 2 fractions in bulk
	frac_1 = c_1/c_tot
	frac_2 = c_2/c_tot
	frac_H = c_H/c_tot

	# Solves for potential and ion distrubutions in diffuse region
	def solver(sigma_d, c_1, c_2, x_end, psi_guess=None):

		# Calculated values
		rho_1_bulk = c_1*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_2_bulk = c_2*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_H_bulk = c_H*1000*N_A 		# bulk number density of H+		(m^-3)
		rho_Cl_bulk = c_tot*1000*N_A 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(e**2*(2*rho_Cl_bulk)/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi[0]/(k*T))
			rho_H = rho_H_bulk*np.exp(-e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_1+z_Li*rho_2+rho_H+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
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
		rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi/(k*T))
		rho_H = rho_H_bulk*np.exp(-e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_1plate(x_dist, -(x_2+x_1), -x_2)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		sigma_0, sigma_beta_1, sigma_beta_2, psi_0 = sigma_0_beta_calc_2cation_pH(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, L, C_1, C_2, psi_beta, eng)
		sigma_beta = sigma_beta_1+sigma_beta_2
		psi = append_1plate(psi, psi_0, psi_beta)
		dpsi = append_1plate(dpsi, 0, 0)
		rho_1 = append_1plate(rho_1, 0, 0)
		rho_2 = append_1plate(rho_2, 0, 0)
		rho_H = append_1plate(rho_H, 0, 0)
		rho_Cl = append_1plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_1, rho_2, rho_H, rho_Cl, dpsi))
		return Solution_2cation_pH(solution, sigma_0, sigma_beta, sigma_beta_1, sigma_beta_2, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c_1, c_2, x_end, psi_guess=None):
		solution = solver(sigma_d, c_1, c_2, x_end, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	# Sets overflow errors to warning
	np.seterr(over='warn')
	warnings.filterwarnings('error')

	# Print progress header
	print_1plate_2cation_ads_pH_header(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, C_1, C_2)

	guess = 0

	# Initialize guess
	guess_exist = False
	while not guess_exist:
		try:
			# guess = sigma_d_guess_helper(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, C_1, C_2, solver_1plate_1cation_ads_pH)
			sigma_d_guess = fsolve(lambda sigma_d: objective(sigma_d, c_1, c_2, x_end), guess)
			solution = solver(sigma_d_guess, c_1, c_2, x_end)
			psi_guess = solution.psi_guess_1plate()
			guess_exist = True
			print_1plate_ads(x_end, sigma_d_guess)
		except Warning:
			x_end = x_end/2

	# Iterates until dpsi/dx = 0 as x -> inf
	while not solution_end(solution):
		x_end_prev = x_end
		x_end += interval
		good = False
		while not good:
			try:
				sigma_d_guess = fsolve(
					lambda sigma_d: objective(sigma_d, c_1, c_2, x_end, psi_guess),
					sigma_d_guess)
				solution = solver(sigma_d_guess, c_1, c_2, x_end, psi_guess)
				print_1plate_ads(x_end, sigma_d_guess)
				psi_guess = solution.psi_guess_1plate()
				good = True
			except Warning:
				x_end = np.mean([x_end_prev, x_end])
	return solution





def solver_2plate_1cation_ads(c, K_ads, D, **kwargs):

	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter_string_list:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	# Solves for potential and ion distrubutions in diffuse region
	def solver(sigma_d, c, psi_guess=False):

		# Calculated values
		rho_Na_bulk = c*1000*N_A 	# bulk number density of Na+	(m^-3)
		rho_Cl_bulk = c*1000*N_A 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(2*e**2*rho_Na_bulk/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_Na+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
		def bc(psia, psib):
			return np.array([psia[1]-sigma_d/(eps_0*eps_bulk), psib[1]+psia[1]])

		size = 50
		x_dist = np.linspace(0, D, size)
		if type(psi_guess) == bool:
			psi_guess = -sigma_d/(eps_0*eps_bulk*kappa)*np.exp(-kappa*x_dist[0:size//2])
			psi_guess = np.concatenate((psi_guess, list(reversed(psi_guess))))
			dpsi_guess = sigma_d/(eps_0*eps_bulk)*np.exp(-kappa*x_dist[0:size//2])
			dpsi_guess = np.concatenate((dpsi_guess, list(reversed(dpsi_guess))))
			psi_guess = np.vstack((psi_guess, dpsi_guess))

		res = solve_bvp(fun, bc, x_dist, psi_guess)
		psi = res.sol(x_dist)[0]
		dpsi = res.sol(x_dist)[1]
		rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_2plate(x_dist, -(x_2+x_1), -x_2, True)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		psi_0 = psi_0_calc(psi_beta, sigma_0, C_1)
		psi = append_2plate(psi, psi_0, psi_beta)
		dpsi = append_2plate(dpsi, 0, 0)
		rho_Na = append_2plate(rho_Na, 0, 0)
		rho_Cl = append_2plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		sigma_beta = sigma_beta_calc_1cation(c, K_ads, L, psi_beta)
		return Solution_1cation(solution, sigma_0, sigma_beta, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c, psi_guess=False):
		solution = solver(sigma_d, c, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	if c >= 1:
		sigma_d = fsolve(lambda sigma_d: objective(sigma_d, c), 0)
		return solver(sigma_d, c)
	else:
		np.seterr(over='warn')
		warnings.filterwarnings('error')

		c_prev = 0
		sigma_d_guess = 0.003
		psi_guess = False
		c_end = 1
		c_half = c*0.5

		print(' ')
		print('-'*85)
		print('Conc (M)\tD (nm)\t\tCurr Conc (M)\t\t\tSigma_d Guess (C/m^2)')
		print('-'*85)

		c_list = []
		sigma_d_list = []
		psi_guess_list = []

		while c_end > 0.6*c:
			c_curr = (c_half+c_end)/2
			good = False
			while not good:
				try:
					sigma_d_guess = fsolve(lambda sigma_d: objective(sigma_d, c_curr, psi_guess), sigma_d_guess)
					solution = solver(sigma_d_guess, c_curr, psi_guess)
					psi_guess = np.vstack((solution.psi[2:-2], solution.dpsi[2:-2]))
					c_list += [c_curr]
					sigma_d_list += [sigma_d_guess[0]]
					psi_guess_list += [psi_guess]
					c_prev = c_end
					c_end = c_curr
					good = True
					print(
						'{0:.1e}\t\t{1:.4}\t\t'.format(c, D*1e9) +
						'{0:.18}\t\t{1}'.format(c_curr, sigma_d_guess))
				except Warning:
					c_curr = (c_curr + c_end)/2

		c_list.reverse()
		sigma_d_list.reverse()
		psi_guess_list.reverse()
		sigma_d_func = interp1d(c_list, sigma_d_list, 'cubic')
		sigma_d = sigma_d_func(c)
		psi_guess = np.zeros(psi_guess_list[0].shape)
		
		m = 0
		for _ in psi_guess:
			n = 0
			for _ in psi_guess[m]:
				entry_list = np.zeros(len(psi_guess_list))
				entry_index = 0
				for entry in psi_guess_list:
					entry_list[entry_index] = entry[m][n]
					entry_index += 1
				entry_func = interp1d(c_list, entry_list)
				psi_guess[m][n] = entry_func(c)
				n += 1
			m += 1

		return solver(sigma_d, c, psi_guess)





def solver_2plate_2cation_ads(c_1, c_2, K_ads_1, K_ads_2, D, **kwargs):

	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter_string_list:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	# Total cation concentration
	c_tot = c_1+c_2

	# Computes cation 1 and cation 2 fractions in bulk
	frac_1 = c_1/c_tot
	frac_2 = c_2/c_tot

	# Solves for potential and ion distrubutions in diffuse region
	def solver(sigma_d, c_1, c_2, psi_guess=False):

		# Calculated values
		rho_1_bulk = c_1*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_2_bulk = c_2*1000*N_A 		# bulk number density of Na+	(m^-3)
		rho_Cl_bulk = c_tot*1000*N_A 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(e**2*(rho_1_bulk+rho_2_bulk+rho_Cl_bulk)/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_1+z_Li*rho_2+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
		def bc(psia, psib):
			return np.array([psia[1]-sigma_d/(eps_0*eps_bulk), psib[1]+psia[1]])

		size = 50
		x_dist = np.linspace(0, D, size)
		if type(psi_guess) == bool:
			psi_guess = -sigma_d/(eps_0*eps_bulk*kappa)*np.exp(-kappa*x_dist[0:size//2])
			psi_guess = np.concatenate((psi_guess, list(reversed(psi_guess))))
			dpsi_guess = sigma_d/(eps_0*eps_bulk)*np.exp(-kappa*x_dist[0:size//2])
			dpsi_guess = np.concatenate((dpsi_guess, list(reversed(dpsi_guess))))
			psi_guess = np.vstack((psi_guess, dpsi_guess))

		res = solve_bvp(fun, bc, x_dist, psi_guess)
		psi = res.sol(x_dist)[0]
		dpsi = res.sol(x_dist)[1]
		rho_1 = rho_1_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_2 = rho_2_bulk*np.exp(-z_Li*e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_2plate(x_dist, -(x_2+x_1), -x_2, True)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		psi_0 = psi_0_calc(psi_beta, sigma_0, C_1)
		psi = append_2plate(psi, psi_0, psi_beta)
		dpsi = append_2plate(dpsi, 0, 0)
		rho_1 = append_2plate(rho_1, 0, 0)
		rho_2 = append_2plate(rho_2, 0, 0)
		rho_Cl = append_2plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_1, rho_2, rho_Cl, dpsi))
		sigma_beta, sigma_beta_1, sigma_beta_2 = sigma_beta_calc_2cation(c_1, c_2, K_ads_1, K_ads_2, L, psi_beta)
		return Solution_2cation(solution, sigma_0, sigma_beta, sigma_beta_1, sigma_beta_2, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c_1, c_2, psi_guess=False):
		solution = solver(sigma_d, c_1, c_2, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	if c_tot >= 2:
		sigma_d = fsolve(lambda sigma_d: objective(sigma_d, c_1, c_2), 0)
		return solver(sigma_d, c_1, c_2)
	else:
		np.seterr(over='warn')
		warnings.filterwarnings('error')

		c_tot_prev = 0
		sigma_d_guess = 0.003
		psi_guess = False
		c_tot_end = 2
		c_tot_half = c_tot*0.5

		print(' ')
		print('-'*93)
		print('Conc 1 (M)\tConc 2 (M)\tD (nm)\tCurr Tot Conc (M)\t\tSigma_d Guess (C/m^2)')
		print('-'*93)

		c_tot_list = []
		sigma_d_list = []
		psi_guess_list = []

		while c_tot_end > 0.6*c_tot:
			c_tot_curr = (c_tot_half+c_tot_end)/2
			good = False
			while not good:
				try:
					c_1_curr = c_tot_curr*frac_1
					c_2_curr = c_tot_curr*frac_2
					sigma_d_guess = fsolve(
						lambda sigma_d: objective(sigma_d, c_1_curr, c_2_curr, psi_guess),
						sigma_d_guess)
					solution = solver(sigma_d_guess, c_1_curr, c_2_curr, psi_guess)
					psi_guess = np.vstack((solution.psi[2:-2], solution.dpsi[2:-2]))
					c_tot_list += [c_tot_curr]
					sigma_d_list += [sigma_d_guess[0]]
					psi_guess_list += [psi_guess]
					c_tot_prev = c_tot_end
					c_tot_end = c_tot_curr
					good = True
					print(
						'{0:.1e}\t\t{1:.1e}\t\t{2:.4}\t'.format(c_1, c_2, D*1e9) +
						'{0:.18}\t\t{1}'.format(c_tot_end, sigma_d_guess))
				except Warning:
					c_tot_curr = (c_tot_curr + c_tot_end)/2

		c_tot_list.reverse()
		sigma_d_list.reverse()
		psi_guess_list.reverse()
		sigma_d_func = interp1d(c_tot_list, sigma_d_list, 'cubic')
		sigma_d = sigma_d_func(c_tot)
		psi_guess = np.zeros(psi_guess_list[0].shape)
		
		m = 0
		for _ in psi_guess:
			n = 0
			for _ in psi_guess[m]:
				entry_list = np.zeros(len(psi_guess_list))
				entry_index = 0
				for entry in psi_guess_list:
					entry_list[entry_index] = entry[m][n]
					entry_index += 1
				entry_func = interp1d(c_tot_list, entry_list)
				psi_guess[m][n] = entry_func(c_tot)
				n += 1
			m += 1

		return solver(sigma_d, c_1, c_2, psi_guess)





def solver_2plate_1cation_ads_pH(c, pH, K_ads, pKa, D, **kwargs):

	if kwargs is not None:
		for key, value in kwargs.items():
			if key in parameter:
				parameter[key] = value

	C_1 = parameter['C_1']
	C_2 = parameter['C_2']

	# Solves for potential and ion distrubutions in diffuse region
	def solver(sigma_d, c, psi_guess=None):

		# Calculated values
		rho_Na_bulk = c*1000*N_A 	# bulk number density of Na+	(m^-3)
		rho_H_bulk = 10**(-pH)*1000*N_A
		rho_Cl_bulk = rho_Na_bulk+rho_H_bulk 	# bulk number density of Cl-	(m^-3)
		kappa = np.sqrt(2*e**2*(rho_Na_bulk+rho_H_bulk+rho_Cl_bulk)/(eps_0*eps_bulk*k*T))

		# Poisson-Boltzmann equation system
		def fun(x, psi):
			rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi[0]/(k*T))
			rho_H = rho_H_bulk*np.exp(-e*psi[0]/(k*T))
			rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi[0]/(k*T))
			dpsi_dx = psi[1]
			d2psi_dx2 = -e/(eps_0*eps_bulk)*(z_Na*rho_Na+rho_H+z_Cl*rho_Cl)
			return np.vstack((dpsi_dx, d2psi_dx2))

		# Boundary conditions
		def bc(psia, psib):
			return np.array([psia[1]-sigma_d/(eps_0*eps_bulk), psib[1]+psia[1]])

		size = 50
		x_dist = np.linspace(0, D, size)
		if psi_guess is None:
			psi_guess = -sigma_d/(eps_0*eps_bulk*kappa)*np.exp(-kappa*x_dist[0:size//2])
			psi_guess = np.concatenate((psi_guess, list(reversed(psi_guess))))
			dpsi_guess = sigma_d/(eps_0*eps_bulk)*np.exp(-kappa*x_dist[0:size//2])
			dpsi_guess = np.concatenate((dpsi_guess, list(reversed(dpsi_guess))))
			psi_guess = np.vstack((psi_guess, dpsi_guess))

		res = solve_bvp(fun, bc, x_dist, psi_guess)
		psi = res.sol(x_dist)[0]
		dpsi = res.sol(x_dist)[1]
		rho_Na = rho_Na_bulk*np.exp(-z_Na*e*psi/(k*T))
		rho_H = rho_H_bulk*np.exp(-e*psi/(k*T))
		rho_Cl = rho_Cl_bulk*np.exp(-z_Cl*e*psi/(k*T))
		x_dist = append_2plate(x_dist, -(x_2+x_1), -x_2, True)
		psi_beta = psi_beta_calc(psi[0], sigma_d, C_2)
		sigma_0, sigma_beta, psi_0 = sigma_0_beta_calc_1cation_pH(c, pH, K_ads, pKa, L, C_1, C_2, psi_beta)
		psi = append_2plate(psi, psi_0, psi_beta)
		dpsi = append_2plate(dpsi, 0, 0)
		rho_Na = append_2plate(rho_Na, 0, 0)
		rho_H = append_2plate(rho_H, 0, 0)
		rho_Cl = append_2plate(rho_Cl, 0, 0)

		# return np.vstack((x_dist, psi, rho_Na, rho_Cl, dpsi))
		solution = np.vstack((x_dist, psi, rho_Na, rho_H, rho_Cl, dpsi))
		# sigma_0, sigma_beta = sigma_0_beta_calc_1cation_pH(c, pH, K_ads, pKa, L, psi_beta)
		return Solution_1cation_pH(solution, sigma_0, sigma_beta, sigma_d)

	# Finds initial condition for solver
	def objective(sigma_d, c, psi_guess=None):
		solution = solver(sigma_d, c, psi_guess)
		sigma_0 = solution.sigma_0
		sigma_beta = solution.sigma_beta
		sigma_d = solution.sigma_d
		return sigma_0+sigma_beta+sigma_d

	if c >= 1:
		sigma_d = fsolve(lambda sigma_d: objective(sigma_d, c), 0)
		return solver(sigma_d, c)
	else:

		# Sets overflow errors to warning
		np.seterr(over='warn')
		warnings.filterwarnings('error')

		# Print header for output display
		print_2plate_1cation_ads_pH_header(c, pH, K_ads, pKa, C_1, C_2, D)

		# Initialize guess at 1 M
		c_curr = 1.0
		sigma_d_guess = fsolve(lambda sigma_d: objective(sigma_d, c_curr), 0)
		solution = solver(sigma_d_guess, 1)
		psi_guess = np.vstack((solution.psi[2:-2], solution.dpsi[2:-2]))

		# Iterate through concentrations until final concentration is reached
		while c_curr > c:
			if c_curr/c - 1 < 0.1:
				try:
					sigma_d = fsolve(
						lambda sigma_d: objective(sigma_d, c, psi_guess),
						sigma_d_guess)
					print_2plate_1cation_ads_pH(c, sigma_d)
					return solver(sigma_d, c, psi_guess)
				except Warning:
					pass
			c_end = c_curr
			c_curr = np.mean((c_curr, c))
			good = False
			while not good:
				try:
					sigma_d_guess = fsolve(
						lambda sigma_d: objective(sigma_d, c_curr, psi_guess),
						sigma_d_guess)
					solution = solver(sigma_d_guess, c_curr, psi_guess)
					psi_guess = np.vstack((solution.psi[2:-2], solution.dpsi[2:-2]))
					good = True
					print_2plate_1cation_ads_pH(c_curr, sigma_d_guess)
				except Warning:
					c_curr = np.mean((c_curr, c_end))





