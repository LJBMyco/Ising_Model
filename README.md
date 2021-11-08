# Ising_Model
Python simulation of the Ising Model on a 2D Lattice
Lattice sites can be selected using either Kawaski or Glauber dynamics
The decision to swap the polarity of said site is then made using the Metropolis alogrithim

Usage: ising.py Method (Kawasaki/Glauber) Lattice Size Temperature Sweeps Data/Animate

Data: 
Total Energy, Total Magnetisation, Heat Capacity and Suseptability can be collected for every tenth sweep, after equilibriating at 100 sweeps.
These are calculated at a temperature range of 1 to 3 in steps of 0.1
ising_plot.py can then be used to create graphs from the output data
Temperature is irrelavent

Animate:
Display the output indefinetly, with each fram representing a full sweep
Sweeps is irrelavent
