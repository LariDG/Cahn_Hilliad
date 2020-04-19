"""
Modelling and Visualisation in Physics
Checkpoint 3: PDEs
Class to initialise 2D potential and electric fields 
of a lattice with a singular monopole at the centre.
Done so based on Poisson statistics.
Author: L. Dorman-Gajic
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import math 
from scipy import signal


class Poisson2D(object):

    def __init__(self, size, omega, thres):
        """
            initialising 

            :param size: size of 2d lattice as a tuple
            :param omega: over relaxation parameter
            :param thres: the threshold for convergence         
        """
        self.size = size
        self.omega = omega
        self.thres = thres

        self.build()

    def build(self):
        """
        Building lattice for the potential field (setting all boundary points to zero)
        with random noise as well as initialising the charge distribution.
        """
        phi_size = (self.size[0]-2, self.size[1]-2)
        self.phi = (np.random.choice(a=[0.01,-0.01], size = phi_size)*np.random.random(phi_size))
        self.phi = np.insert(self.phi,phi_size[0]-2,0,axis=0)
        self.phi = np.insert(self.phi,phi_size[1]-2,0,axis=1)
        self.phi = np.insert(self.phi,0,0,axis=0)
        self.phi = np.insert(self.phi,0,0,axis=1)

        self.rho = np.zeros(self.size)

    def monopole(self):
        """
        Setting up a monopole in the centre of the lattice of the charge distribution.
        """
        self.rho[self.size[0]//2, self.size[1]//2] = 1.0

    def gauss_seidel(self):
        """
        Gauss-Seidel algorithm using for loops for a 2D lattice
        using over relaxation to increase convergence time.
        """
        for i in range(1,self.size[0]-1):
            for j in range(1,self.size[1]-1):
                self.phi[(i,j)] = ((1/4)*(self.phi[(i+1,j)] + self.phi[(i-1,j)] + self.phi[(i,j+1)] + self.phi[(i,j-1)] + self.rho[(i,j)]) - self.phi[(i,j)])*self.omega + self.phi[(i,j)]


    def terminate_condition(self, p_a, p_b):
        """
        if statment to determine if the lattice has converged. 
        """
        if np.sum(abs(p_a - p_b), axis = None) <= self.thres:
            return True
        else:
            return False

        



    
