#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains plotting procedures for plotting the thermodynamic
quantities around a phase transition
"""

import numpy as np
import matplotlib.pyplot as plt

# Set fontsizes in figures
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)

L = np.array([40, 60, 100, 160])
t = np.loadtxt("./Results/evolution_L=40.txt", usecols=0)

E = []
M = []
Cv = []
X = []

for l in L:
    file = "./Results/evolution_L=%s.txt" % (l)
    E.append(np.loadtxt(file, usecols=1))
    M.append(np.loadtxt(file, usecols=2))
    Cv.append(np.loadtxt(file, usecols=3))
    X.append(np.loadtxt(file, usecols=4))


maxX = np.argmax(X, axis=1)
fit = np.polyfit(1. / L, t[maxX], 1)
polynom = np.poly1d(fit)

fig = plt.figure()
plt.plot(1. / L, t[maxX], "o")
plt.plot(1. / L, polynom(1. / L))
plt.xlabel("$1/L$")
plt.ylabel("$T_c$ $[\,J/k_B\,]$")
plt.legend(["Finite Lattice", "%.4f x + %.4f" % (fit[0], fit[1])])
plt.grid()
fig.savefig("./Plots/critTemp.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, E[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\langle E \\rangle /L^2$ $[\,J\,]$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("./Plots/evolution_energy.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, M[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\langle |M| \\rangle /L^2$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("./Plots/evolution_magnetization.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, Cv[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$ C_V /L^2$ $[\,k_B\,]$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("./Plots/evolution_cv.pdf")

fig = plt.figure()
for i in range(len(L)):
    plt.plot(t, X[i], linewidth=0.8)
plt.xlabel("T $[\,J/k_B\,]$")
plt.ylabel("$\\chi /L^2$")
plt.legend(["40x40 Lattice", "60x60 Lattice", "100x100 Lattice",
            "160x160 Lattice"])
plt.grid()
fig.savefig("./Plots/evolution_susceptibility.pdf")
