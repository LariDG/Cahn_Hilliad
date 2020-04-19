# MVP Checkpoint3: Partial Differential Equation Solver

Author: Larisa Dorman-Gajic

        s1752550 
        
        University of Edinburgh 
        
        School of Physics and Astronomy
        
A series of classes:
- Cahn-Hilliard 
- Poisson 
- Poisson_Mag
- Poisson2D
And a main module which runs these classes:
- partial_diff.py

These modules us Cahn-Hilliard algorithm and Poisson statistics to calculate mixtures of oil and water as well as field around 
points and wires.



To run command line should read as:

```python3 partial_diff.py {size of square lattice} {type of simulation to run} {initial conditions for phi: either 0.0, 0.5, -0.5}```

Example dats and plots for all combinations as shown in files "data" and "plots"


## partial_diff.py

The main module which uses classes to run different algorithms based on user input.
This includes an animation of an oil and water mixture based on starting parameters
using the Cahn-Hilliard algorthim as well as being able to produce a graph for the 
free energy density based on this mixture and initial parameters. Also using Poisson
statistics to test the potential, electric and magnetic fields for lattices with 
either a monopole or wire in the centre. Also produces a plot for the convergences of
the potential field of a monopole in a 2d lattice with changing amounts of over-relaxation.

## Cahn-Hilliard

Class using the Cahn-Hilliard algorithm to model an oil and water mixture
using a 2D lattice and animated using Func Animation.
Also calculates the free energy density of the system.

## Poisson

Class to initialise 3D potential and electric fields 
of a lattice with a singular monopole at the centre.
Done so based on Poisson statistics.

## Poisson_Mag

Class to initialise 3D potential and magnetic fields 
of a lattice with a wire through the centre.
Done so based on Poisson statistics.

## Poisson2D

Class to initialise 2D potential and electric fields 
of a lattice with a singular monopole at the centre.
Done so based on Poisson statistics.
