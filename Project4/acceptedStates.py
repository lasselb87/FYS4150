#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains plotting procedures for plotting the number of
accepted states
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Set fontsizes in figures
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)

accepted = []
N = [1000, 2000, 3000, 4000, 5000, 6000]

for n in N:
    os.system("mpirun -np 1 ./simulation.x %s 0 20 2 1" % n)
    accepted.append(np.loadtxt("./Results/meta.txt", usecols=0)[5])

fig = plt.figure()
plt.plot(N, accepted, "o-")
plt.xlabel("Monte Carlo Cycles [#]")
plt.ylabel("Accepted States")
plt.legend(["20x20 Lattice\nT = 2 $[\,J/k_B\,]$"])
plt.gcf().set_tight_layout(True)
plt.grid()
fig.savefig("./Plots/acceptedCycles.pdf")

accepted = []
T = np.linspace(1, 3, 20)
for t in T:
    os.system("mpirun -np 1 ./simulation.x 10000 0 20 %s 1" % t)
    accepted.append(np.loadtxt("./Results/meta.txt", usecols=0)[5])

accepted = np.array(accepted) / 10000

fig = plt.figure()
plt.plot(T, accepted, "o-")
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("Accepted States per Sweep")
plt.legend(["20x20 Lattice\n 10000 cycles"])
plt.gcf().set_tight_layout(True)
plt.grid()
fig.savefig("./Plots/acceptedTemp.pdf")
