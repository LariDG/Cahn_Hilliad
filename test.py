import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from partial_diff import Cahn_Halliard
import sys
import math 


def main():

    size_input = int(sys.argv[1])
    simulate = sys.argv[2]
    phi_0 = float(sys.argv[3])
    size = (size_input, size_input)

    if simulate == "Y":  
        lattice = Cahn_Halliard(size, 1.0, 1.0, 0.1, 0.1, 0.1, phi_0)
        lattice.run(1000000, 1)

    if simulate == "N":
        file_handle = open("free_energy.dat", "w+")
        lattice = Cahn_Halliard(size, 1.0, 1.0, 0.1, 0.1, 0.1, phi_0)
        sweeps = 100000

        time = []
        fed = []

        for i in range(sweeps):
            if i % 50 == 0:
                fed.append(np.sum(lattice.free_energy()))
                print(i)
                time.append(i)

            lattice.update_phi()
            lattice.update_mu()

        for j in range(len(time)):
            file_handle.write(str(time[j]) + " " + str(fed[j]) + "\n")
        file_handle.close()


main()
