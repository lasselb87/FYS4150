#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains procedures for running the simulation and writing the
results to file for multiple temperatures around the critical temperature
"""

import os
import numpy as np
import math as m
import sys

cycles = int(sys.argv[1])
cutoff = int(sys.argv[2])
L = int(sys.argv[3])
T_start = float(sys.argv[4])
T_end = float(sys.argv[5])
T_step = float(sys.argv[6])

N = int(m.ceil((T_end - T_start) / T_step))

file = open("./Results/evolution_L=%s.txt" % L, "w")

for i in range(N + 1):
    os.system("mpirun -np 2 ./simulation.x %s %s %s %s %s" %
              (cycles, cutoff, L, (T_start + i * T_step), i))
    expect = np.loadtxt("./Results/expectation.txt", usecols=0)
    for val in expect:
        file.write("%s " % val)
    file.write("\n")
    print("go!")

file.close()
