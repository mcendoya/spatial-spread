# -*- coding: utf-8 -*-
"""

Initial introduction at t=0 for the different types of initial introduction 
(random, 5 foci and one focus)
100 inital infected individuals
"""

import functions_spread as fn
import pandas as pd
import numpy as np
from random import sample


# Data with coordinates of individuals
xy = pd.read_csv("pts_v2.csv")
# Coordinates
coords = xy[['x', 'y']].values


# RANDOM

for x in range(1, 11):
    R = np.zeros(len(xy))
    idx = sample(range(0, len(R)), 100)
    R[idx] = 1
    xy["R_{0}".format(x)] = R


# 5 FOCI

# focus radious (m)
D_MAX = 1000

for x in range(1, 11):
    f5_t0 = np.zeros(len(xy))
    idx = sample(range(0, len(f5_t0)), 5)
    f5_t0[idx] = 1
    indI = idx
    indS = np.where(f5_t0 == 0)[0]
    nb_1000 = fn.nb_dist(D_MAX, coords, indI, indS)
    nb_I = []
    for i in range(len(nb_1000)):
        idxI = sample(nb_1000[i], 19)
        nb_I.append(idxI)
    nb_I_flat = [j for i in nb_I for j in i]
    f5_t0[nb_I_flat] = 1
    xy["f5_{0}".format(x)] = f5_t0


# 1 FOCUS

# focus radious (m)
D_MAX = 1000

for x in range(1, 11):
    f1_t0 = np.zeros(len(xy))
    idx = sample(range(0, len(f1_t0)), 1)
    f1_t0[idx] = 1
    indI = idx
    indS = np.where(f1_t0 == 0)[0]
    nb_1000 = fn.nb_dist(D_MAX, coords, indI, indS)
    nb_I = sample(nb_1000[0], 99)
    f1_t0[nb_I] = 1
    xy["f1_{0}".format(x)] = f1_t0

# SAVE INTRODUCTION DATA
xy.to_csv('coords.csv', index=False)
