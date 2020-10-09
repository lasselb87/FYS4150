#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# Set fontsizes in figures
params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large'}
plt.rcParams.update(params)


def u(x):
    '''Closed form solution'''
    return 1 - (1 - np.exp(-10)) * x - np.exp(-10 * x)


def plotResult(n, x_num, v, title):
    x_exact = np.linspace(0, 1, 10000)
    plt.plot(x_exact, u(x_exact), 'k', lw=3.5, label='Analytic')
    plt.plot(x_num, v, 'r', lw=1.3, label='Numerical with n=' +
             str(n) + ' gridpoints')
    plt.grid()
    plt.gca().set_xlabel('$x$')
    plt.gca().set_ylabel('$u(x)$')
    plt.gca().set_title(title)
    plt.legend(loc='best')
    plt.show()


def plot_error(h_err, eps_err, title):
    plt.plot(np.log10(h_err), eps_err, label='Relative Error')
    plt.plot(np.log10(h_err), 2 * np.log10(h_err),
             label='Reference: Graph with slope=2')
    plt.grid()
    plt.gca().set_xlabel('$\\log(h)$')
    plt.gca().set_ylabel('$\epsilon$')
    plt.gca().set_title(title)
    plt.legend(loc='best')
    plt.show()


file1 = '/Users/Nicolai/Documents/Atom/CompPhys/Project1/project1_results.txt'
x_num = np.loadtxt(file1, usecols=0)
v1 = np.loadtxt(file1, usecols=1)       # General Algorithm result
v2 = np.loadtxt(file1, usecols=2)       # Optimized Algorithm result
v3 = np.loadtxt(file1, usecols=3)       # LU-decomp Algorithm result

file2 = '/Users/Nicolai/Documents/Atom/CompPhys/Project1/RelativeError.txt'
h_err = np.loadtxt(file2, usecols=0)
eps_err = np.loadtxt(file2, usecols=1)

n = 1000      # n = 10, 100, 1000
# Plot result from General Algorithm
plotResult(n, x_num, v1, 'General Algorithm')
# Plot result from Optimized Algorithm
plotResult(n, x_num, v2, 'Optimized Algorithm')
# Plot result from LU-decomp Algorithm
plotResult(n, x_num, v3, 'LU-decomposition Algorithm')
# Plot relative error
plot_error(h_err, eps_err,
           'log-log Plot of Relative Error as a \n Function of Step-size')
