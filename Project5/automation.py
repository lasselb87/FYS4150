#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program contains procedures for compiling and running all the programs in
this project, reproducing all results
"""

import os

# Compile programs and run tests
os.system("make test.x")
os.system("make main.x")
os.system("./main.x")
os.system("python3 plot.py")
