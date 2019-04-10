"""
SOLUTION CLASSES
================

This module contains several classes used to turn solver output into an object.

"""

import numpy as np

class Solution_1cation:
	"""This is used for systems with 1 type of cation.

	"""
	def __init__(self, solution, sigma_0, sigma_beta, sigma_d):
		self.x = solution[0]
		self.psi = solution[1]
		self.cation = solution[2]
		self.anion = solution[3]
		self.dpsi = solution[4]
		self.sigma_beta = sigma_beta
		self.sigma_d = float(sigma_d)
		self.sigma_0 = sigma_0

class Solution_2cation:
	"""This is used for systems with 2 types of cations.

	"""
	def __init__(self, solution, sigma_0, sigma_beta, sigma_beta_1, sigma_beta_2, sigma_d):
		self.x = solution[0]
		self.psi = solution[1]
		self.cation_1 = solution[2]
		self.cation_2 = solution[3]
		self.anion = solution[4]
		self.dpsi = solution[5]
		self.sigma_beta = sigma_beta
		self.sigma_beta_1 = sigma_beta_1
		self.sigma_beta_2 = sigma_beta_2
		self.sigma_d = float(sigma_d)
		self.sigma_0 = sigma_0

	def psi_guess_1plate(self):
		return np.vstack((self.psi[2:], self.dpsi[2:]))

	def psi_guess_2plate(self):
		return np.vstack((self.psi[2:-2], self.dpsi[2:-2]))

class Solution_1cation_pH:
	"""This is used for systems with 1 type of cation with pH effects.

	"""
	def __init__(self, solution, sigma_0, sigma_beta, sigma_d):
		self.x = solution[0]
		self.psi = solution[1]
		self.cation = solution[2]
		self.H = solution[3]
		self.anion = solution[4]
		self.dpsi = solution[5]
		self.sigma_beta = sigma_beta
		self.sigma_d = float(sigma_d)
		self.sigma_0 = sigma_0

	def psi_guess_1plate(self):
		return np.vstack((self.psi[2:], self.dpsi[2:]))

	def psi_guess_2plate(self):
		return np.vstack((self.psi[2:-2], self.dpsi[2:-2]))

class Solution_2cation_pH:
	"""This is used for systems with 2 type of cation with pH effects.

	"""
	def __init__(self, solution, sigma_0, sigma_beta, sigma_beta_1, sigma_beta_2, sigma_d):
		self.x = solution[0]
		self.psi = solution[1]
		self.cation_1 = solution[2]
		self.cation_2 = solution[3]
		self.H = solution[4]
		self.anion = solution[5]
		self.dpsi = solution[6]
		self.sigma_beta = sigma_beta
		self.sigma_beta_1 = sigma_beta_1
		self.sigma_beta_2 = sigma_beta_2
		self.sigma_d = float(sigma_d)
		self.sigma_0 = sigma_0

	def psi_guess_1plate(self):
		return np.vstack((self.psi[2:], self.dpsi[2:]))

	def psi_guess_2plate(self):
		return np.vstack((self.psi[2:-2], self.dpsi[2:-2]))