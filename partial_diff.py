import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
import math 

class Cahn_Halliard(object):
    def __init__(self, size, dx, dt, a, M, kappa, phi_0):
        self.size = size
        self.dx = dx
        self.dt = dt
        self.a = a
        self.M = M
        self.kappa = kappa
        self.phi_0 = phi_0

        self.build(phi_0)


    def build(self, phi_0):
        self.phi = (np.random.choice(a=[1.0,-1.0], size = self.size)*np.random.random(self.size)) + phi_0
        #self.phi = np.zeros(self.size)
        #for i in range(self.size[0]):
         #   for j in range(self.size[1]):
          #      self.phi[i, j] = np.random.uniform(-0.1, 0.11) + self.phi_0
        self.mu = np.zeros(self.size)

    def pbc(self, indices):
        return(indices[0]%self.size[0], indices[1]%self.size[1])  

    #def mu(self, indices):
      #  mu = (-self.a * self.phi[indices] + self.a * self.phi[indices]**3 - self.kappa * self.laplacian(self.phi, indices))
       # return mu

    def laplacian(self, lattice):
        lap = np.zeros(self.size)
        for i in range(self.size[0]):
            for j in range(self.size[1]):
                #indices = (i,j)
                laplacian_x = (lattice[self.pbc((i+1,j))] + lattice[self.pbc((i-1,j))] - 2 * lattice[i,j]) / self.dx**2
                laplacian_y = (lattice[self.pbc((i,j+1))] + lattice[self.pbc((i,j-1))] - 2 * lattice[i,j]) / self.dx**2
                lap[i,j] = laplacian_x + laplacian_y
                #lap[i,j] = (lattice[self.pbc((i+1,j))] + lattice[self.pbc((i-1,j))] + lattice[self.pbc((i,j+1))] + lattice[self.pbc((i,j-1))] - (4 * lattice[self.pbc((i,j))]))/self.dx**2
        return lap
    
    def differential(self, lattice):
        diff = np.zeros(self.size)
        for i in range(self.size[0]):
            for j in range(self.size[1]):
                diff[i,j] = (lattice[self.pbc((i+1,j))] - lattice[self.pbc((i-1,j))] + lattice[self.pbc((i,j+1))] - lattice[self.pbc((i,j-1))])/self.dx*2.0
        print("diff")
        print(diff)
        return diff

    def update_phi(self):
        #phi = np.zeros(self.size)
        #for n in range(self.size[0]):
            #for m in range(self.size[1]):
                #indices = (n,m)
        lap = self.laplacian(self.mu)
        self.phi = self.phi + (self.M * self.dt) * lap
                

    def update_mu(self):
        #for n in range(self.size[0]):
         #   for m in range(self.size[1]):
        lap = self.laplacian(self.phi)
        self.mu = (-self.a * self.phi + self.a * self.phi**3 - self.kappa * lap)

    def free_energy(self):

        #diff = self.differential(self.phi)
        diff = np.gradient(self.phi)
        print(self.phi)
        print("fe")
        free_energy = -(self.a/2.0) * self.phi**2 + (self.a/4.0) * self.phi**4 + (self.kappa/2.0) * (diff[0]**2 + diff[1]**2)
        print(free_energy)
        return free_energy
       
    def run(self, iterations, it_per_frame):
        """
        method running the data into FuncAnimation
        """
        self.it_per_frame = it_per_frame
        self.figure = plt.figure()
        self.image = plt.imshow(self.phi, cmap='seismic', animated=True, interpolation='gaussian')
        self.animation = animation.FuncAnimation(self.figure, self.animate, repeat=False, frames=iterations, interval=20 , blit=True)
        plt.clim(-1.1, 1.1)
        plt.colorbar()
        plt.show()

    def animate(self, *args):
        """

        a loop to but data into animation
        """
        for t in range(100):
            self.update_phi()
            self.update_mu()
        self.image.set_array(self.phi)
        return self.image,