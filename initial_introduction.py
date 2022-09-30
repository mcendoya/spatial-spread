# -*- coding: utf-8 -*-
"""
@author: Martina Cendoya

Initial introduction at t=0 in the different scenarios: 
random and aggregated (one or five foci), 
in the whole area and in the test area.
"""

from random import sample
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import model_functions as fn

#######################################################################
#########################   WHOLE AREA   ##############################
#######################################################################

XY = pd.read_csv("./Data/coords_alm.csv")
COORDS = XY[['x', 'y']].values

#################################################

# 100 random initial infected


R_t0 = np.zeros(len(XY))
idx = sample(range(0, len(R_t0)), 100)
R_t0[idx] = 1


# 100 initial infected in 1 focus

# focus radious (m)
D_MAX = 1000

f1_t0 = np.zeros(len(XY))
idx = sample(range(0, len(f1_t0)), 1)
f1_t0[idx] = 1
indI = idx
indS = np.where(f1_t0 == 0)[0]
nb_1000 = fn.nb_dist(D_MAX, COORDS, indI, indS)
nb_I = sample(nb_1000[0], 99)
f1_t0[nb_I] = 1

# 100 initial infected in 5 foci

f5_t0 = np.zeros(len(XY))
idx = sample(range(0, len(f5_t0)), 5)
f5_t0[idx] = 1
indI = idx
indS = np.where(f5_t0 == 0)[0]
nb_1000 = fn.nb_dist(D_MAX, COORDS, indI, indS)
nb_I = []
for i in range(len(nb_1000)):
    idxI = sample(nb_1000[i], 19)
    nb_I.append(idxI)
nb_I_flat = [j for i in nb_I for j in i]
f5_t0[nb_I_flat] = 1

xy = XY[['x', 'y']]
xy['R_t0'] = R_t0
xy['f1_t0'] = f1_t0
xy['f5_t0'] = f5_t0


plt.scatter(xy['x'], xy['y'], color="green", s=10)
plt.scatter(xy['x'][xy['R_t0'] == 1], xy['y'][xy['R_t0'] == 1], color="red", s=10)
plt.show()

plt.scatter(xy['x'], xy['y'], color="green", s=10)
plt.scatter(xy['x'][xy['f1_t0'] == 1], xy['y'][xy['f1_t0'] == 1], color="red", s=10)
plt.show()

plt.scatter(xy['x'], xy['y'], color="green", s=10)
plt.scatter(xy['x'][xy['f5_t0'] == 1], xy['y'][xy['f5_t0'] == 1], color="red", s=10)
plt.show()

xy.to_csv('./Data/coords.csv', index=False)

#######################################################################
##########################   TEST AREA   ##############################
#######################################################################

# 25 cells
cells = list(range(599, 603+1))
cells.extend(list(range(565, 569+1)))
cells.extend(list(range(532, 536+1)))
cells.extend(list(range(500, 504+1)))
cells.extend(list(range(459, 463+1)))
XY_EX = XY[XY['cell'].isin(cells)]

COORDS = XY_EX[['x', 'y']].values

# 10 random initial infected

R_t0 = np.zeros(len(XY_EX))
idx = sample(range(0, len(R_t0)), 10)
R_t0[idx] = 1


# 10 initial infected in 1 focus

# focus radious (m)
D_MAX = 1000

f1_t0 = np.zeros(len(XY_EX))
idx = sample(range(0, len(f1_t0)), 1)
f1_t0[idx] = 1
indI = idx
indS = np.where(f1_t0 == 0)[0]
nb_1000 = fn.nb_dist(D_MAX, COORDS, indI, indS)
nb_I = sample(nb_1000[0], 9)
f1_t0[nb_I] = 1

xy_ex = XY_EX[['x', 'y']]
xy_ex['R_t0'] = R_t0
xy_ex['f1_t0'] = f1_t0

plt.scatter(xy_ex['x'], xy_ex['y'], color="green", s=10)
plt.scatter(xy_ex['x'][xy_ex['R_t0'] == 1], xy_ex['y'][xy_ex['R_t0'] == 1], color="red", s=10)
plt.show()

plt.scatter(xy_ex['x'], xy_ex['y'], color="green", s=10)
plt.scatter(xy_ex['x'][xy_ex['f1_t0'] == 1], xy_ex['y'][xy_ex['f1_t0'] == 1], color="red", s=10)
plt.show()


xy_ex.to_csv('./Data/coords_ex.csv', index=False)

###################################################
