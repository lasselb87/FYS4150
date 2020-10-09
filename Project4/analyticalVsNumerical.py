#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains plotting procedures for plotting numerical results
compared with analytical solutions for the 2x2 lattice Ising model
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# Set fontsizes in figures
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)

N = sys.argv[1]

file = "./Results/evolution_L=2.txt"
t = np.loadtxt(file, usecols=0)
T = np.linspace(t[0], t[-1], 1000)
E_anal = -(8 * np.sinh(8 / T)) / (np.cosh(8 / T) + 3) / 4
M_anal = (2 * np.exp(8 / T) + 4) / (np.cosh(8 / T) + 3) / 4
Cv_anal = 1 / T**2 * 64 * (3 * np.cosh(8 / T) + 1) / \
    (np.cosh(8 / T) + 3)**2 / 4
X_anal = 1 / T * ((8 * np.exp(8 / T) + 8) /
                  (np.cosh(8 / T) + 3) - (4 * M_anal)**2) / 4


E = np.loadtxt(file, usecols=1)
M = np.loadtxt(file, usecols=2)
Cv = np.loadtxt(file, usecols=3)
X = np.loadtxt(file, usecols=4)

fig = plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.xlabel("$T$ $[\,J/k_B\,]$")
plt.ylabel("$\\langle E \\rangle /L^2$ $[\,J\,]$")
plt.plot(T, E_anal, linewidth=1.5)
plt.plot(t, E)
plt.legend(["Analytical", "Numerical"], loc="best")
plt.grid()

plt.subplot(2, 2, 2)
plt.xlabel("$T$ $[\,J/k_B\,]$")
plt.ylabel("$\\langle |M| \\rangle /L^2$")
plt.plot(T, M_anal, linewidth=1.5)
plt.plot(t, M)
plt.legend(["Analytical", "Numerical"], loc="best")
plt.grid()

plt.subplot(2, 2, 3)
plt.xlabel("$T$ $[\,J/k_B\,]$")
plt.ylabel("$ C_V /L^2$ $[\,k_B\,]$")
plt.plot(T, Cv_anal, linewidth=1.5)
plt.plot(t, Cv)
plt.legend(["Analytical", "Numerical"], loc="best")
plt.grid()

plt.subplot(2, 2, 4)
plt.xlabel("$T$ $[\,J/k_B\,]$")
plt.ylabel("$\\chi /L^2$")
plt.plot(T, X_anal, linewidth=1.5)
plt.plot(t, X)
plt.legend(["Analytical", "Numerical"], loc="best")
plt.grid()

plt.gcf().set_tight_layout(True)
fig.savefig("./Plots/numericalVsAnalytical_N=%s.pdf" % N)
