"""
Modelling and Visualisation in Physics
Checkpoint 3: PDEs
The main module which uses classes to run different algorithms based on user input.
This includes an animation of an oil and water mixture based on starting parameters
using the Cahn-Hilliard algorthim as well as being able to produce a graph for the 
free energy density based on this mixture and initial parameters. Also using Poisson
statistics to test the potential, electric and magnetic fields for lattices with 
either a monopole or wire in the centre. Also produces a plot for the convergences of
the potential field of a monopole in a 2d lattice with changing amounts of over-relaxation.
Author: L. Dorman-Gajic
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from partial_diff import Cahn_Halliard
from Poisson import Poisson
from Poisson_Mag import Poisson_Mag
from Poisson2D import Poisson2D
import sys
import math 


def main():
    #input by user of the size of square lattice and the type algorithm they wish to produce
    size_input = int(sys.argv[1])
    simulate = sys.argv[2]

    if simulate == "Y":  
        #This simulates the Cahn-Hilliard algorithm
        #user input of initial condition of mixture 
        phi_0 = float(sys.argv[3])
        size = (size_input, size_input)
        #calling class to initialise lattice
        lattice = Cahn_Halliard(size, 1.0, 1.0, 0.1, 0.1, 0.1, phi_0)
        #runs simulation in class
        lattice.run(1000000, 1)

    elif simulate == "CH":
        #to produce free energy density graph
        #user input of initial condition of mixture 
        phi_0 = float(sys.argv[3])
        size = (size_input, size_input)
        #open file
        file_handle = open("free_energy.dat", "w+")
        #initialise lattice using class
        lattice = Cahn_Halliard(size, 1.0, 1.0, 0.1, 0.1, 0.1, phi_0)
        sweeps = 1000000

        #staring lists to plot 
        time = []
        fed = []

        for i in range(sweeps):
            #only taking every 100th sweep
            if i % 100 == 0:
                #calculating free energy density
                fed.append(np.sum(lattice.free_energy()))
                print(i)
                time.append(i)
            #updating the lattices each step
            lattice.update_phi()
            lattice.update_mu()

        #writing to file
        for j in range(len(time)):
            file_handle.write(str(time[j]) + " " + str(fed[j]) + "\n")
        file_handle.close()

    elif simulate == "Poisson":
        #producing lattices for potential and electric field, 3D
        thres = 0.01
        size = (size_input, size_input, size_input)
        #initialising lattice in class
        lattice = Poisson(size, thres)
        sweeps = 1000000

        #creating monopole
        lattice.monopole()

        #condition for ending the run
        end = False


        for i in range(sweeps):
            #storing initial phi value for determining convergence
            before = np.array(lattice.phi)
            #updating phi
            lattice.phi = lattice.jacobi(lattice.phi)
            if lattice.terminate_condition(before, lattice.phi) == False:
                print(i)
                pass
            elif lattice.terminate_condition(before, lattice.phi) == True:
                #end when convergence occurs
                print("Jacobi converged in {} steps".format(i))
                break

        #storing field in data and text files for plotting
        phi_xy_plane = lattice.phi[:,:,25]
        np.savetxt("Potential_Field.txt", phi_xy_plane)
        with open("Potential_Field.dat", "w+") as dataFile:

            for i in range(size[0]):
                for j in range(size[1]):
                    dataFile.write('%lf, %lf, %lf\n' %(float(i), float(j), phi_xy_plane[(i,j)]))
                

        #storing E-field in data files for plotting
        E_x, E_y, E_z = lattice.e_field()
        E_x_plane = E_x[:,:,25]
        E_y_plane = E_y[:,:,25]


        with open("Electric_Field.dat", "w+") as dataFile:
            for i in range(size[0]):
                for j in range(size[1]):
                    dataFile.write('%lf, %lf, %lf, %lf\n' %(float(i), float(j), E_x_plane[(i,j)], E_y_plane[(i,j)]))



    elif simulate == "Poisson_Mag":
        #producing lattices for potential and magnetic field, 3D
        thres = 0.01
        size = (size_input, size_input, size_input)
        #initialising lattice in class
        lattice = Poisson_Mag(size, 0.0, thres)
        sweeps = 1000000

        #creating wire
        lattice.wire()

        end = False


        for i in range(sweeps):
            #storing initial value for convergence condition
            before = np.array(lattice.A)
            #updating lattice
            lattice.A = lattice.jacobi(lattice.A)
            if lattice.terminate_condition(before, lattice.A) == False:
                print(i)
                pass
            elif lattice.terminate_condition(before, lattice.A) == True:
                #end when convergence occurs
                print("Jacobi converged in {} steps".format(i))
                break

        #plotting potential field
        A_xy_plane = lattice.A[:,:,25]
        np.savetxt("Potential_Field_Mag.txt", A_xy_plane)
        with open("Potential_Field_Mag.dat", "w+") as dataFile:

            for i in range(size[0]):
                for j in range(size[1]):
                    dataFile.write('%lf, %lf, %lf\n' %(float(i), float(j), A_xy_plane[(i,j)]))
                

        #plotting magnetic field
        B_x, B_y, B_z = lattice.m_field()
        B_x_plane = B_x[:,:,25]
        B_y_plane = B_y[:,:,25]


        with open("Magnetic_Field.dat", "w+") as dataFile:
            for i in range(size[0]):
                for j in range(size[1]):
                    dataFile.write('%lf, %lf, %lf, %lf\n' %(float(i), float(j), B_x_plane[(i,j)], B_y_plane[(i,j)]))
 

    elif simulate == "SOR":
        #over relaxation for 2d poisson lattice
        phi_0 = float(sys.argv[3])
        thres = 0.01
        size = (size_input, size_input)
        sweeps = 1000000
        file_handle = open("omega_sweeps.dat", "w+")

        #array of omega (over relaxation) values
        omega_values = np.arange(1.0, 2.0, 0.01)
        #initialising list for number of sweeps till convergence occurs
        converge_sweep = []
        for i in range(len(omega_values)):
            omega = omega_values[i]
            #initialise lattice from class using varrying omega values
            lattice = Poisson2D(size, 1.0, 1.0, 0.0, omega, thres)
            lattice.monopole()

            end = False

            for j in range(sweeps):
                #storing initial lattice so to determine convergence 
                before = np.array(lattice.phi)
                lattice.gauss_seidel()
                if lattice.terminate_condition(before, lattice.phi) == False:
                    print(j)
                    pass
                elif lattice.terminate_condition(before, lattice.phi) == True:
                    #end and store sweep if convergence occured
                    converge_sweep.append(j)
                    break 

        #save data file for sweeps and omegas to plot 
        for m in range(len(converge_sweep)):
            file_handle.write(str(omega_values[m]) + " " + str(converge_sweep[m]) + "\n")
        file_handle.close()
        

main()