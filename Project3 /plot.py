#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains plotting procedures for producing the results in
this project
"""

import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats

# Set fontsizes in figures
params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
plt.rcParams.update(params)


def curvefit(x, y, deg=1):
    coefs = poly.polyfit(x, y, deg)
    ffit = poly.Polynomial(coefs)    # instead of np.poly1d
    slope = coefs[1]
    return ffit, slope


# -------------------------------------------------------------------------
# Plot the orbit of Earth in the Earth-Sun system
# -------------------------------------------------------------------------
if sys.argv[1] == "singlePlanet":
    file = "./Raw_Data/data.txt"

    name = sys.argv[2]
    if name == "Euler":
        labelname = "Orbit with Forward Euler"
    elif name == "Verlet":
        labelname = "Orbit with Velocity Verlet"

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure(figsize=(10, 10))
    plt.gca().set_aspect("equal")
    plt.plot([0, 0], [0, 0], 'o', color='orange', label='Sun', markersize=25)
    plt.plot(x[0], y[0], 'o', color='deepskyblue',
             label='Earth, initial position', mfc='none', markersize=10)
    plt.plot(x, y, ':', color='deepskyblue', label=labelname)
    plt.plot(x[-1], y[-1], 'o', color='deepskyblue',
             label='Earth, final position', markersize=10)
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.legend(loc='upper right')
    fig.savefig("./Results/earthOrbit_%s.pdf" % (name))

# -------------------------------------------------------------------------
# Plot benchmark of CPU times
# -------------------------------------------------------------------------
if sys.argv[1] == "benchmark":
    file = "./Raw_Data/benchmark.txt"
    N = np.log10(np.loadtxt(file, usecols=0))
    T_euler = np.log10(np.loadtxt(file, usecols=1))
    T_verlet = np.log10(np.loadtxt(file, usecols=2))
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(
        N, T_euler)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(
        N, T_verlet)

    fig = plt.figure(figsize=(10, 5))
    plt.plot(N, T_euler, 'o:', label="Forward Euler method")
    plt.plot(N, intercept1 + slope1 * N, '-',
             label='Euler fitted line.\nSlope:' + str(round(slope1, 4)))
    plt.plot(N, T_verlet, 'o:', label="Velocity Verlet method")
    plt.plot(N, intercept2 + slope2 * N, '-',
             label='Verlet fitted line.\nSlope: ' + str(round(slope2, 4)))
    plt.gca().set_xlabel('$\\log_{10}(N)$')
    plt.gca().set_ylabel('$\\log_{10}(T_{CPU})$')
    plt.subplots_adjust(right=0.68)
    plt.legend(loc='center left', bbox_to_anchor=(1.03, 0.5),
               fancybox=True, borderaxespad=0, ncol=1)
    plt.grid(True)
    fig.savefig("./Results/benchmark.pdf")

# -------------------------------------------------------------------------
# Plot total energy as a function of time
# -------------------------------------------------------------------------
if sys.argv[1] == "energy":
    file = "./Raw_Data/energy.txt"

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure(figsize=(10, 5))
    plt.plot(np.linspace(0, 5, len(x)), x, label='Forward Euler method')
    plt.plot(np.linspace(0, 5, len(y)), y, label='Velocity Verlet method')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.gca().set_xlabel('$T$ [yr]')
    plt.gca().set_ylabel('$E_{TOT}$ [$M_\\odot$AU$^2$/yr$^2$]')
    plt.grid(True)
    plt.legend(loc='best')
    fig.savefig("./Results/energyFluctuation.pdf")

# -------------------------------------------------------------------------
# Plot fluctuation as a function of the number of intergration points
# -------------------------------------------------------------------------
if sys.argv[1] == "fluctuation":
    file = "./Raw_Data/fluctuation.txt"
    label_ = sys.argv[2]
    if label_ == "Euler":
        labelname = "Fluctation with Forward Euler"
    elif label_ == "Verlet":
        labelname = "Fluctuation with Velocity Verlet"
    elif label_ == "EarthAndJupiter":
        labelname = "Fluctuation with Velocity Verlet"

    n = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure()
    plt.plot(np.log10(n), np.log10(y), label=labelname)
    plt.gca().set_xlabel('$\\log_{10}(N)$')
    plt.gca().set_ylabel('$\\log_{10}(\\epsilon)$')
    plt.grid(True)
    plt.legend(loc='best')
    fig.savefig("./Results/fluctuation_%s.pdf" % label_)

    if label_ == "Euler" or label_ == "Verlet":
        y2 = np.loadtxt(file, usecols=2)
        fig = plt.figure()
        plt.plot(np.log10(n), np.log10(y2))
        plt.gca().set_xlabel('$\\log_{10}(N)$')
        plt.gca().set_ylabel('$\\log_{10}(\\eta)$')
        plt.grid(True)
        fig.savefig("results/fluctuation_angular_%s.pdf" % label_)


# -------------------------------------------------------------------------
# Plot different initial velocities to find escape velocity
# -------------------------------------------------------------------------
if sys.argv[1] == "escapeVel":
    mylegend = []
    fig = plt.figure(figsize=(10, 5))
    for i in range(6):
        file = "./Raw_Data/data" + str(i) + ".txt"
        x = np.loadtxt(file, usecols=0)
        y = np.loadtxt(file, usecols=1)
        labelname = "$v_0$ = $\\sqrt{%s}$$2\\pi$" % (1 + 0.25 * i)
        plt.plot([0, 0], [0, 0], 'o', color='orange', markersize=15)
        plt.plot(x, y, label=labelname)
    plt.xlabel("$X$ [AU]")
    plt.ylabel("$Y$ [AU]")
    plt.grid(True)
    # plt.legend(mylegend)
    plt.legend(loc='best')
    fig.savefig("./Results/escapeVel.pdf")


# -------------------------------------------------------------------------
# Plot inverse-square -> inverse-cube gravity
# -------------------------------------------------------------------------
if sys.argv[1] == "beta":
    file = "./Raw_Data/data.txt"

    name = sys.argv[2]
    if name == "beta1":
        labelname = "$1/r^{2.5}$"
    elif name == "beta2":
        labelname = "$1/r^3$"

    x = np.loadtxt(file, usecols=0)
    y = np.loadtxt(file, usecols=1)

    fig = plt.figure(figsize=(10, 10))
    plt.gca().set_aspect("equal")
    plt.plot([0, 0], [0, 0], 'o', color='orange', label='Sun', markersize=25)
    plt.plot(x[0], y[0], 'o', color='deepskyblue',
             label='Earth, initial position', mfc='none', markersize=10)
    plt.plot(x, y, ':', color='deepskyblue', label=labelname)
    plt.plot(x[-1], y[-1], 'o', color='deepskyblue',
             label='Earth, final position', markersize=10)
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.legend(loc='upper right')
    fig.savefig("./Results/earthBeta_%s.pdf" % (name))


# -------------------------------------------------------------------------
# Plot the motion of the Earth-Jupiter-Sun system
# -------------------------------------------------------------------------
if sys.argv[1] == "earthAndJupiter":
    for i in range(1, 4):
        fig = plt.figure(figsize=(10, 10))
        x1 = np.loadtxt("./Raw_Data/data%s.txt" % i, usecols=0)
        y1 = np.loadtxt("./Raw_Data/data%s.txt" % i, usecols=1)

        x2 = np.loadtxt("./Raw_Data/data%s.txt" % i, usecols=3)
        y2 = np.loadtxt("./Raw_Data/data%s.txt" % i, usecols=4)

        plt.plot([0, 0], [0, 0], 'o', color='orange',
                 label='Sun', markersize=27)

        plt.plot(x1[0], y1[0], 'o', color='deepskyblue',
                 label='Earth, initial position', mfc='none', markersize=6)
        plt.plot(x1, y1, ':', color='deepskyblue', label="Earth's orbit")
        plt.plot(x1[-1], y1[-1], 'o', color='deepskyblue',
                 label='Earth, final position', markersize=6)
        plt.plot(x2[0], y2[0], 'o', color='chocolate',
                 label='Jupiter, initial position', mfc='none',
                 markersize=9 * i)
        plt.plot(x2, y2, ':', color='chocolate', label="Jupiter's orbit")
        plt.plot(x2[-1], y2[-1], 'o', color='chocolate',
                 label='Jupiter, final position', markersize=9 * i)
        plt.xlabel("$X$ [AU]")
        plt.ylabel("$Y$ [AU]")
        plt.legend(loc='best')
        fig.savefig("./Results/earthAndJupiter%s.pdf" % i)


# -------------------------------------------------------------------------
# Plot the motion of all planets including the sun
# -------------------------------------------------------------------------
if sys.argv[1] == "all":
    file = "./Raw_Data/data.txt"

    pos_sun = np.loadtxt(file, usecols=(0, 1, 2))
    pos_mercury = np.loadtxt(file, usecols=(3, 4, 5))
    pos_venus = np.loadtxt(file, usecols=(6, 7, 8))
    pos_earth = np.loadtxt(file, usecols=(9, 10, 11))
    pos_mars = np.loadtxt(file, usecols=(12, 13, 14))
    pos_jupiter = np.loadtxt(file, usecols=(15, 16, 17))
    pos_saturn = np.loadtxt(file, usecols=(18, 19, 20))
    pos_uranus = np.loadtxt(file, usecols=(21, 22, 23))
    pos_neptune = np.loadtxt(file, usecols=(24, 25, 26))
    pos_pluto = np.loadtxt(file, usecols=(27, 28, 29))

    # All planets
    fig = plt.figure(figsize=(15, 15))
    plt.gcf().add_subplot(111, projection='3d')
    plt.gca().set_aspect("equal")
    plt.grid(False)
    plt.plot(pos_sun[:, 0], pos_sun[:, 1],
             pos_sun[:, 2], 'yo', label='Sun', markersize=5)
    plt.plot(pos_mercury[:, 0], pos_mercury[:, 1],
             pos_mercury[:, 2], color='darkslategrey', label='Mercury')
    plt.plot(pos_venus[:, 0], pos_venus[:, 1],
             pos_venus[:, 2], color='tan', label='Venus')
    plt.plot(pos_earth[:, 0], pos_earth[:, 1],
             pos_earth[:, 2], color='deepskyblue', label='Earth')
    plt.plot(pos_mars[:, 0], pos_mars[:, 1],
             pos_mars[:, 2], color='tomato', label='Mars')
    plt.plot(pos_jupiter[:, 0], pos_jupiter[:, 1],
             pos_jupiter[:, 2], color='orange', label='Jupiter')
    plt.plot(pos_saturn[:, 0], pos_saturn[:, 1],
             pos_saturn[:, 2], color='olive', label='Saturn')
    plt.plot(pos_uranus[:, 0], pos_uranus[:, 1],
             pos_uranus[:, 2], color='darkviolet', label='Uranus')
    plt.plot(pos_neptune[:, 0], pos_neptune[:, 1],
             pos_neptune[:, 2], color='darkblue', label='Neptune')
    plt.plot(pos_pluto[:, 0], pos_pluto[:, 1],
             pos_pluto[:, 2], color='black', label='Pluto')
    plt.gca().dist = 7.2
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.gca().set_zlabel('$Z$ [AU]')
    plt.legend(loc='best', prop={'size': 25})
    fig.savefig("./Results/allPlanets.png")

    # Inner planets
    fig = plt.figure(figsize=(10, 10))
    plt.gcf().add_subplot(111, projection='3d')
    plt.gca().set_aspect("equal")
    plt.grid(False)
    plt.plot(pos_sun[:, 0], pos_sun[:, 1],
             pos_sun[:, 2], 'yo', label='Sun', markersize=5)
    plt.plot(pos_mercury[:, 0], pos_mercury[:, 1],
             pos_mercury[:, 2], color='darkslategrey', label='Mercury')
    plt.plot(pos_venus[:, 0], pos_venus[:, 1],
             pos_venus[:, 2], color='tan', label='Venus')
    plt.plot(pos_earth[:, 0], pos_earth[:, 1],
             pos_earth[:, 2], color='deepskyblue', label='Earth')
    plt.plot(pos_mars[:, 0], pos_mars[:, 1],
             pos_mars[:, 2], color='tomato', label='Mars')
    plt.gca().dist = 7.7
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.gca().set_zlabel('$Z$ [AU]')
    plt.gca().set_zlim(-1.5, 1.5)
    plt.legend(loc='best', prop={'size': 20})
    fig.savefig("./Results/innerPlanets.png")

    # Outer planets
    fig = plt.figure(figsize=(15, 15))
    plt.gcf().add_subplot(111, projection='3d')
    plt.gca().set_aspect("equal")
    plt.grid(False)
    plt.plot(pos_sun[:, 0], pos_sun[:, 1],
             pos_sun[:, 2], 'yo', label='Sun', markersize=15)
    plt.plot(pos_jupiter[:, 0], pos_jupiter[:, 1],
             pos_jupiter[:, 2], color='orange', label='Jupiter')
    plt.plot(pos_saturn[:, 0], pos_saturn[:, 1],
             pos_saturn[:, 2], color='olive', label='Saturn')
    plt.plot(pos_uranus[:, 0], pos_uranus[:, 1],
             pos_uranus[:, 2], color='darkviolet', label='Uranus')
    plt.plot(pos_neptune[:, 0], pos_neptune[:, 1],
             pos_neptune[:, 2], color='darkblue', label='Neptune')
    plt.gca().dist = 7.7
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.gca().set_zlabel('$Z$ [AU]')
    plt.gca().set_zlim(-2.5, 2.5)
    plt.legend(loc='best', prop={'size': 20})
    fig.savefig("./Results/outerPlanets.png")

    # Earth-Jupiter-Sun
    fig = plt.figure(figsize=(10, 10))
    plt.gcf().add_subplot(111, projection='3d')
    plt.gca().set_aspect("equal")
    plt.grid(False)
    plt.plot(pos_sun[:, 0], pos_sun[:, 1],
             pos_sun[:, 2], 'yo', label='Sun', markersize=15)
    plt.plot(pos_earth[:, 0], pos_earth[:, 1],
             pos_earth[:, 2], color='deepskyblue', label='Earth')
    plt.plot(pos_jupiter[:, 0], pos_jupiter[:, 1],
             pos_jupiter[:, 2], color='orange', label='Jupiter')
    plt.gca().dist = 7.7
    plt.gca().set_xlabel('$X$ [AU]')
    plt.gca().set_ylabel('$Y$ [AU]')
    plt.gca().set_zlabel('$Z$ [AU]')
    plt.gca().set_zlim(-2.5, 2.5)
    plt.legend(loc='best', prop={'size': 20})
    fig.savefig("./Results/threeBodyDynamic.png")


# -------------------------------------------------------------------------
# Plot the precession of the perihelion of Mercury
# -------------------------------------------------------------------------
if sys.argv[1] == "precession":

    rad_to_arcsec = (360 * 60 * 60) / (2 * np.pi)
    file1 = "./Raw_Data/angle1.txt"
    radians = np.arctan(np.loadtxt(file1, usecols=0))
    arcseconds = radians * rad_to_arcsec
    t = np.linspace(0, 100, len(arcseconds))
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(
        t, arcseconds)

    file2 = "./Raw_Data/angle2.txt"
    radians2 = np.arctan(np.loadtxt(file2, usecols=0))
    arcseconds2 = radians2 * rad_to_arcsec
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(
        t, arcseconds2)

    fig = plt.figure(figsize=(10, 5))
    plt.plot(t, arcseconds2, label="Einsten precession")
    plt.plot(t, intercept2 + slope2 * t, 'k', lw=2,
             label='Einstein fitted line.\nSlope:' + str(round(slope2, 4)))
    plt.plot(t, arcseconds, label="Newton precession")
    plt.plot(t, intercept1 + slope1 * t, 'r', lw=2,
             label='Newton fitted line.\nSlope:' + str(round(slope1, 4)))
    plt.gca().set_xlabel("$T$ [yr]")
    plt.gca().set_ylabel("$\\theta_p$ [arcseconds]")
    plt.legend(loc='best')
    fig.savefig("./Results/precession.pdf")
