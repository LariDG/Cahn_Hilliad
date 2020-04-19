"""
Modelling and Visualisation in Physics
Checkpoint 3: PDEs
Class to initialise 3D potential and magnetic fields 
of a lattice with a wire through the centre.
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


class Poisson_Mag(object):

    def __init__(self, size, A_0, thres):
        """
            initialising 

            :param size: size of 3d lattice as a tuple
            :param A_0: initial condition
            :param thres: the threshold for convergence         
        """
        self.size = size
        self.omega = 1.0
        self.A_0 = A_0
        self.thres = thres

        self.build()

    def build(self):
        """
        building lattices for A (the vector potential) and J (the current field)
        """
        A_size = (self.size[0]-2, self.size[1]-2, self.size[2]-2)
        self.A = (np.random.choice(a=[0.01,-0.01], size = A_size)*np.random.random(A_size) + self.A_0)
        self.A = np.insert(self.A,A_size[0]-2,0,axis=0)
        self.A = np.insert(self.A,A_size[1]-2,0,axis=1)
        self.A = np.insert(self.A,A_size[2]-2,0,axis=2)
        self.A = np.insert(self.A,0,0,axis=0)
        self.A = np.insert(self.A,0,0,axis=1)
        self.A = np.insert(self.A,0,0,axis=2)

        self.J = np.zeros(self.size)

    def wire(self):
        """
        putting a wire through the centre of the current field
        """
        self.J[self.size[0]//2, self.size[1]//2, :] = 1.0 / self.size[2]

    def jacobi(self, lattice):
        """
        Convolution method of updating via the jacobi algorithm given a lattice 
        """
        kernel = np.array([[[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,0.0]],
                           [[0.0,1.0,0.0],[1.0,0.0,1.0],[0.0,1.0,0.0]],
                           [[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,0.0]]])
        return ((signal.fftconvolve(lattice, kernel, mode='same')  + self.J)/ 6.0)

    def gauss_seidel(self):
        """
        Gauss-Seidel algorithm using for loops for a 3D lattice.
        """
        for i in range(1,self.size[0]-1):
            for j in range(1,self.size[1]-1):
                for k in range(1,self.size[2]-1):
                    self.A[(i,j,k)] = ((1/6)*(self.A[(i+1,j,k)] + self.A[(i-1,j,k)] + self.A[(i,j+1,k)] + self.A[(i,j-1,k)] + self.A[(i,j,k+1)] + self.A[(i,j,k-1)] + self.J[(i,j,k)]) - self.A[(i,j,k)])*self.omega + self.A_0[(i,j,k)]

    def m_field(self):
        """
        Calculating the magnetic field from the gradient of the potential field
        """
        grad = np.gradient(self.A)

        B_x = grad[1] - grad[2]
        B_y = - grad[2] - grad[0]
        B_z = - grad[0] - grad[1]
        return (B_x, B_y, B_z)
    
        
    def terminate_condition(self, p_a, p_b):
        """
        if statment to determine if the lattice has converged. 
        """
        if np.sum(abs(p_a - p_b), axis = None) <= self.thres:
            return True
        else:
            return False

