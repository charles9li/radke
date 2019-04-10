"""
HELPER
======

This module contains a list of helper functions.

"""

import os, sys
from scipy.integrate import simps
from scipy.optimize import root
from solvers import *

##########
# APPEND #
##########

def append_1plate(list_old, item_1, item_2):
	"""Returns a new numpy array that has the numpy array [item_1, item_2]
	appended to the beginning of the old list. Used in 1-plate solvers.

	Examples
	--------
	>>> import numpy as np
	>>> append_1plate([3, 4], 1, 2)
	array([1, 2, 3, 4])

	"""
	return np.append(np.array([item_1, item_2]), list_old)

def append_2plate(list_old, item_1, item_2, is_x=False):
	"""Returns a new numpy array that has the numpy arrays [item_1, item_2] and
	[item_2, item_1] appended to the beginning and end of the old list,
	respectively. However, if the list is for the x-coordinate
	(i.e., is_x=True), then the numpy array [-item_2, -item_1]+list_old[-1] is
	appended to the end of the old list instead. Used in 2-plate solvers.

	Examples
	--------
	>>> import numpy as np
	>>> append_2plate([3, 3], 1, 2)
	array([1, 2, 3, 3, 2, 1])
	>>> append_2plate([0, 1, 2], -2, -1, True)
	array([-2, -1,  0,  1,  2,  3,  4])

	"""
	list_new = np.append(np.array([item_1, item_2]), list_old)
	if is_x:
		return np.append(list_new, np.array([-item_2, -item_1])+list_new[-1])
	else:
		return np.append(list_new, np.array([item_2, item_1]))


#######################
# 1 PLATE END CHECKER #
#######################

def solution_end(solution):
	"""This function checks to see if the termination condition is satisfied for
	1-plate solvers.

	"""
	psi_2 = solution.psi[-2]
	psi_1 = solution.psi[-1]
	return abs(psi_2-psi_1) < 1e-5


#############
# POTENTIAL #
#############

def psi_0_calc(psi_beta, sigma_0, C_1):
	"""This function computes 0-plane potential given beta-plane potential,
	surface charge density, and integral capacitance between the 0-plane and
	the beta-plane.

	"""
	return psi_beta + sigma_0/C_1

def psi_beta_calc(psi_d, sigma_d, C_2):
	"""This function computes beta-plane potential given d-plane potential,
	diffuse layer charge density, and integral capacitance between the beta-
	plane and the d-plane.

	"""
	return psi_d - sigma_d/C_2


##################
# CHARGE DENSITY #
##################

def sigma_beta_calc_1cation(c_bulk, K_ads, L, psi_beta):
	"""This function computes beta-plane charge density given bulk cation
	concentration, K_ads, surface charge density, beta-plane potential. Used
	for 1-cation systems.

	"""
	return e*L*K_ads*c_bulk/(K_ads*c_bulk+np.exp(e*psi_beta/(k*T)))

def sigma_beta_calc_2cation(c_1, c_2, K_ads_1, K_ads_2, L, psi_beta):
	"""This function computes beta-plane charge density given bulk cation
	concentration, K_ads, surface charge density, beta-plane potential. Used
	for 2-cation systems.

	"""
	psi_red = np.exp(e*psi_beta/(k*T))
	A = K_ads_1*c_1+psi_red
	B = K_ads_2*c_2+psi_red
	SA = K_ads_1*c_1*L*psi_red/(A*B-K_ads_1*K_ads_2*c_1*c_2)
	SB = K_ads_2*c_2/B*(L-SA)
	return e*(SA+SB), e*SA, e*SB

def sigma_0_beta_calc_1cation_pH(c_bulk, pH, K_ads, pKa, L, C_1, C_2, psi_beta):
	"""This function computes beta-plane charge density given bulk cation
	concentration, K_ads, pKa, pH, surface charge density, beta-plane potential. 
	Used for 1-cation systems.

	"""
	psi_beta_red = float(np.exp(e*psi_beta/(k*T)))
	c_H = 10**(-pH)
	K_a = 10**(-float(pKa))

	def equations(p):
		sigma_0, sigma_beta, psi_0 = p
		psi_diff_red = psi_beta_red/np.exp(e*psi_0/(k*T))
		SM = K_ads*c_bulk*L/(K_ads*c_bulk+c_H/K_a*psi_diff_red+psi_beta_red)
		S = SM*psi_beta_red/(K_ads*c_bulk)
		return (psi_0-psi_beta-sigma_0/C_1, sigma_0+e*(SM+S), sigma_beta-e*SM)

	return fsolve(equations, (1, 1, 1))

test = False
if test:
	def sigma_0_beta_calc_2cation_pH(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, L, C_1, C_2, psi_beta, eng):
		"""This function computes beta-plane charge density given bulk cation
		concentration, K_ads, pKa, pH, surface charge density, beta-plane potential. 
		Used for 1-cation systems.
		
		"""
		c_1 = float(c_1)
		c_1 = eng.double(c_1)
		c_2 = float(c_2)
		c_2 = eng.double(c_2)
		K_ads_1 = float(K_ads_1)
		K_ads_1 = eng.double(K_ads_1)
		K_ads_2 = float(K_ads_2)
		K_ads_2 = eng.double(K_ads_2)
		pKa = float(pKa)
		pKa = eng.double(pKa)
		psi_beta = float(psi_beta[0])
		psi_beta = eng.double(psi_beta)
		S, SM1, SM2, SH, psi_0 = eng.sigma_0_beta_calc_2cation_pH(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, L, C_1, C_2, psi_beta, T, nargout=5)
		S = float(S)
		SM1 = float(SM1)
		SM2 = float(SM2)
		SH = float(SH)
		psi_0 = float(psi_0)
		sigma_0 = -e*(L-SH)
		sigma_beta_1 = e*SM1
		sigma_beta_2 = e*SM2
		return sigma_0, sigma_beta_1, sigma_beta_2, psi_0
else:
	def sigma_0_beta_calc_2cation_pH(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, L, C_1, C_2, psi_beta, eng):
		"""This function computes beta-plane charge density given bulk cation
		concentration, K_ads, pKa, pH, surface charge density, beta-plane potential. 
		Used for 1-cation systems.
		
		"""
		psi_beta_red = float(np.exp(e*psi_beta/(k*T)))
		c_H = 10**(-pH)
		K_a = 10**(-float(pKa))

		def equations(p):
			S, SM_1, SM_2, SH = p
			# S = L/(1 + (K_ads_1*c_1+K_ads_2*c_2)*psi_beta_red + c_H/K_a*psi_0_red)
			sigma_0 = -e*(S+SM_1+SM_2)
			psi_0 = psi_beta + sigma_0/C_2
			return (
				(L-S-SM_1-SM_2-SH)/L,
				(K_ads_1-SM_1/(S*c_1)*psi_beta_red)/K_ads_1,
				(K_ads_2-SM_2/(S*c_2)*psi_beta_red)/K_ads_2,
				(K_a-S*c_H/SH*np.exp(-e*psi_0/(k*T)))/K_a)

		guess = (L*.2, L*.2, L*.2, L*.2)
		# guess = (2.118661739599816e+16, 1.149980827435163e+18, 1.036952955306729e+17, 7.251372596381661e+17)
		# sol = newton_krylov(equations, guess, method='gmres')
		sol = root(equations, guess)
		S, SM_1, SM_2, SH = sol.x
		# print(sol.x)
		# print(sol.success)
		# print(sol.message)
		# scipy root
		# - hybrid/lv
		sigma_0 = -e*(L-SH)
		sigma_beta_1 = e*SM_1
		sigma_beta_2 = e*SM_2
		psi_0 = psi_beta + sigma_0/C_2
		return sigma_0, sigma_beta_1, sigma_beta_2, psi_0


#######################
# SURFACE ION DESNITY #
#######################

def diffuse_ion_conc_2plate(x, rho, rho_0):
	"""
	This function computes the excess surface concentration of an ion located in
	the diffuse region. It takes in a 2-plate profile and computes the density
	for only one side.

	"""
	x = x[0:len(x)//2]
	rho = rho[0:len(rho)//2]
	rho_excess = rho - rho_0
	return simps(rho_excess, x)


##########
# FORCES #
##########

# def pressure_osmotic_1cation():


#########
# PRINT #
#########

def print_1plate_1cation_ads_pH_header(c, pH, K_ads, pKa, C_1, C_2):
	"""
	Prints the progress header for a 1-plate 1-cation system with pH effects.

	"""
	print(' ')
	print('='*70)
	print('Conc\t= %.1e M\t|| C_1 = %.0f microF/cm2' % (c, C_1*100))
	print('pH\t= %.2f\t\t|| C_2 = %.0f microF/cm2' % (pH, C_2*100))
	print('K_ads\t= %.1f M^-1\t||' % (K_ads))
	print('pKa\t= %.2f\t\t||' % (pKa))
	print('='*70)
	print('x_end (nm)\tSigma_d Guess (C/m^2)')
	print('-'*70)

def print_1plate_2cation_ads_header(c_1, c_2, K_ads_1, K_ads_2, C_1, C_2):
	"""
	Prints the progress header for a 1-plate 2-cation system without pH effects.

	"""
	print(' ')
	print('='*70)
	print('conc_1\t= %.1e M\t|| C_1 = %.0f microF/cm2' % (c_1, C_1*100))
	print('conc_2\t= %.1e M\t|| C_2 = %.0f microF/cm2' % (c_2, C_2*100))
	print('K_ads_1\t= %.1f M^-1\t||' % (K_ads_1))
	print('K_ads_2\t= %.1f M^-1\t||' % (K_ads_2))
	print('='*70)
	print('x_end (nm)\tSigma_d Guess (C/m^2)')
	print('-'*70)

def print_1plate_2cation_ads_pH_header(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, C_1, C_2):
	"""
	Prints the progress header for a 1-plate 2-cation system with pH effects.

	"""
	print(' ')
	print('='*70)
	print('conc_1\t= %.1e M\t|| C_1 = %.0f microF/cm2' % (c_1, C_1*100))
	print('conc_2\t= %.1e M\t|| C_2 = %.0f microF/cm2' % (c_2, C_2*100))
	print('pH\t= %.2f\t\t||' % (pH))
	print('K_ads_1\t= %.1f M^-1\t||' % (K_ads_1))
	print('K_ads_2\t= %.1f M^-1\t||' % (K_ads_2))
	print('pKa\t= %.2f\t\t||' % (pKa))
	print('='*70)
	print('x_end (nm)\tSigma_d Guess (C/m^2)')
	print('-'*70)

def print_1plate_ads(x_end, sigma_d_guess):
	"""
	Prints the current cocentration and guess for sigma_d.

	"""
	print('{0:.1f}\t\t{1}'.format(x_end*1e9, sigma_d_guess))

def print_2plate_1cation_ads_pH_header(c, pH, K_ads, pKa, C_1, C_2, D):
	"""
	Prints the progress header for a 2-plate 1-cation system with pH effects.

	"""
	print(' ')
	print('='*70)
	print('Conc\t= %.1e M\t|| C_1 = %.0f microF/cm2' % (c, C_1*100))
	print('pH\t= %.2f\t\t|| C_2 = %.0f microF/cm2' % (pH, C_2*100))
	print('K_ads\t= %.1f M^-1\t||' % (K_ads))
	print('pKa\t= %.2f\t\t||' % (pKa))
	print('D\t= %.1f nm\t||' % (D*1e9))
	print('='*70)
	print('Curr Conc (M)\t\t\tSigma_d Guess (C/m^2)')
	print('-'*70)

def print_2plate_1cation_ads_pH(c_curr, sigma_d_guess):
	"""
	Prints the current cocentration and guess for sigma_d.

	"""
	print('{0:.18}\t\t{1}'.format(c_curr, sigma_d_guess))


#########
# GUESS #
#########

class HiddenPrints:
	def __enter__(self):
		self._original_stdout = sys.stdout
		sys.studout = open(os.devnull, 'w')

	def __exit__(self):
		sys.stdout.close()
		sys.stdout = self._original_stdout

def sigma_d_guess_helper(c_1, c_2, pH, K_ads_1, K_ads_2, pKa, C_1, C_2, solver):
	"""
	This function provides an initial sigma_d guess for the 2-cation pH solvers.
	When the initial guess for sigma_d is too off, it could cause the objective
	function to loop indefinitely.

	"""
	c = c_1 + c_2
	K_ads = (K_ads_1*c_1+K_ads_2+c_2)/c
	solution = solver(c, pH, K_ads, pKa, C_1=C_1, C_2=C_2)
	return solution.sigma_d


###########
# DOCTEST #
###########

if __name__ == "__main__":
	import doctest
	doctest.testmod()