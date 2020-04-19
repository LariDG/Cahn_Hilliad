"""
Modelling and Visualisation in Physics
Checkpoint 3: PDEs
Class to initialise 3D potential and electric fields 
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


class Poisson(object):

    def __init__(self, size, thres):
        """
            initialising 

            :param size: size of 3d lattice as a tuple
            :param thres: the threshold for convergence         
        """
        self.size = size
        self.omega = 1.0
        self.thres = thres

        self.build()


    def build(self):
        """
        Building lattice for the potential field (setting all boundary points to zero)
        with random noise as well as initialising the charge distribution.
        """
        phi_size = (self.size[0]-2, self.size[1]-2, self.size[2]-2)
        self.phi = (np.random.choice(a=[1.0,-1.0], size = phi_size)*np.random.random(phi_size))
        self.phi = np.insert(self.phi,phi_size[0]-2,0,axis=0)
        self.phi = np.insert(self.phi,phi_size[1]-2,0,axis=1)
        self.phi = np.insert(self.phi,phi_size[2]-2,0,axis=2)
        self.phi = np.insert(self.phi,0,0,axis=0)
        self.phi = np.insert(self.phi,0,0,axis=1)
        self.phi = np.insert(self.phi,0,0,axis=2)

        self.rho = np.zeros(self.size)

    def monopole(self):
        """
        Setting up a monopole in the centre of the lattice of the charge distribution.
        """
        self.rho[self.size[0]//2, self.size[1]//2, self.size[2]//2] = 1.0

    def jacobi(self, lattice):
        """
        Convolution method of updating via the jacobi algorithm given a lattice 
        """
        kernel = np.array([[[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,0.0]],
                           [[0.0,1.0,0.0],[1.0,0.0,1.0],[0.0,1.0,0.0]],
                           [[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,0.0]]])
        return ((signal.fftconvolve(lattice, kernel, mode='same')  + self.rho)/ 6.0)


    def gauss_seidel(self):
        """
        Gauss-Seidel algorithm using for loops for a 3D lattice.
        """
        for i in range(1,self.size[0]-1):
            for j in range(1,self.size[1]-1):
                for k in range(1,self.size[2]-1):
                    self.phi[(i,j,k)] = ((1/6)*(self.phi[(i+1,j,k)] + self.phi[(i-1,j,k)] + self.phi[(i,j+1,k)] + self.phi[(i,j-1,k)] + self.phi[(i,j,k+1)] + self.phi[(i,j,k-1)] + self.rho[(i,j,k)]) - self.phi[(i,j,k)])*self.omega + self.phi[(i,j,k)]

    def e_field(self):
        """
        Calculation of the electric field. Returning each component separatly
        """
        E = np.gradient(self.phi)
        return (-E[0], -E[1], -E[2])

    def terminate_condition(self, p_a, p_b):
        """
        if statment to determine if the lattice has converged. 
        """
        if np.sum(abs(p_a - p_b), axis = None) <= self.thres:
            return True
        else:
            return False

    



    

    

    




    

        

