import numpy as np
from scipy.constants import e, k, epsilon_0, N_A


class PashleySolution:

    def __init__(self, c, psi_d, D_list, epsilon=79, T=298, A=2.2e-20):
        self.c = c
        self.psi_d = psi_d
        self.D_list = D_list
        self.epsilon = epsilon*epsilon_0
        self.T = T
        self.A = A
        self.kappa = None
        self.FR_list = None

    def compute_kappa(self):
        rho = 1000*N_A*self.c
        self.kappa = np.sqrt(2*rho*e**2/(self.epsilon*k*self.T))

    def compute_FR(self):
        self.FR_list = np.ones(len(self.D_list))
        if self.kappa is None:
            self.compute_kappa()
        Z = 64*np.pi*self.epsilon*(k*self.T/e)**2*np.tanh(e*self.psi_d/(4*k*self.T))**2
        i = 0
        for D in self.D_list:
            W = (self.kappa**2/(2*np.pi))*Z*np.exp(-self.kappa*D)
            FR = 2*np.pi*W - self.A/(6*D**2)
            self.FR_list[i] = FR
            i += 1
