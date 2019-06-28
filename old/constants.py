"""
CONSTANTS
=========

Provides a number of fundamental constants and parameters specific to modeling
double layer behavior in aqueous systems.

"""


#########################
# PARAMETER STRING LIST #
#########################

parameter_string_list = ['x_1', 'x_2', 'C_1', 'C_2', 'L', 'sigma_0', 'eps_1', 'eps_2']


#########################
# FUNDAMENTAL CONSTANTS #
#########################

k = 1.38e-23		# Boltzmann's constant 			(J/K)
e = 1.6e-19 		# elementary charge 			(C)
N_A = 6.022e23 		# Avogadro's number 			(mol^-1)
eps_0 = 8.85e-12	# permitivitty of free space	(C^2/J*m)
R = 8.314			# gas constant 					(J/mol*K)
F = 96485			# Faraday's constant 			(C/mol)


###################################
# VALENCES FOR CATIONS AND ANIONS #
###################################

z_Li = 1	# valence of Li+
z_Na = 1 	# valence of Na+
z_K = 1		# valence of K+
z_Rb = 1	# valence of Rb+
z_Cs = 1	# valence of Cs+
z_Cl = -1	# valence of Cl-


###############
# IONIC RADII #
###############

R_Li = 60e-12	# Li+ radius	(m)
R_Na = 98e-12	# Na+ radius	(m)
R_K = 133e-12	# K+ radius		(m)
R_Rb = 148e-12	# Rb+ radius	(m)
R_Be = 31e-12	# Be^2+ radius	(m)
R_Mg = 65e-12	# Mg^2+ radius	(m)
R_Ca = 99e-12	# Ca^2+ radius	(m)
R_Sr = 113e-12	# Sr^2+ radius	(m)
R_F = 136e-12	# F- radius		(m)
R_Cl = 181e-12	# Cl- radius	(m)
R_Br = 185e-12	# Br- radius	(m)
R_I = 216e-12	# I- radius		(m)


########################
# IONIC RADII (STOKES) #
########################

R_Na_stokes = 2.76e-10	# Na+ Stokes readius	(m)


##############################
# SYSTEM SPECIFIC PARAMETERS #
##############################

T = 298					# temperature 									(K)
L = 2/(1e-9)**2			# surface site density of mica 					(m^-2)
sigma_0 = -e*L			# surface charge density 						(C/m^2)
eps_1 = 6 				# rel. perm. between wall and IHP
eps_2 = 30				# rel. perm. between IHP and OHP
eps_bulk = 79 			# rel. perm. of bulk (water)
x_1 = R_Na 				# distance from surf to beta plane				(m)
x_2 = R_Na+R_Na_stokes	# distance from beta plane to d plane			(m)
C_1 = eps_1*eps_0/x_1	# capcacitance between surf and beta plane		(F/m^2)
C_2 = eps_2*eps_0/x_2	# capcacitance between beta plane to d plane	(F/m^2)


########################
# PARAMETER DICTIONARY #
########################
parameter = {}
parameter['T'] = 298				# temperature 									(K)
parameter['L'] = 2/(1e-9)**2		# surface site density of mica 					(m^-2)
parameter['sigma_0'] = -e*L			# surface charge density 						(C/m^2)
parameter['C_1'] = eps_1*eps_0/x_1	# capcacitance between surf and beta plane		(F/m^2)
parameter['C_2'] = eps_2*eps_0/x_2	# capcacitance between beta plane to d plane	(F/m^2)