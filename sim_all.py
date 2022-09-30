# -*- coding: utf-8 -*-
"""
@author: Martina Cendoya

Simulations in the whole area
"""

import numpy as np
import pandas as pd
from numba import set_num_threads
import model_functions as fn

set_num_threads(4)

#######################################################################
#########################   WHOLE AREA   ##############################
#######################################################################

XY = pd.read_csv("coords.csv")
COORDS = XY[['x', 'y']].values

RG_LIST = [100, 200, 500, 1000]

np.random.seed(10)
N_NEW = np.random.choice(np.array(range(0, 6)), size=20, replace=True)

#######################################################################

# 1) 100 random initial infected

S_T0 = np.array(XY['R_t0'])

# a) No new events

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015)
    np.savetxt("./results/all_points/R_r"+str(rg)+"_0.txt",
               np.transpose(s1), fmt='%i')

# b) Removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14)
    print(rg)
    np.savetxt("./results/all_points/R_r"+str(rg)+"_r.txt",
               np.transpose(s1), fmt='%i')

# c) New introductions

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/R_r"+str(rg)+"_i.txt",
               np.transpose(s1), fmt='%i')

# d) New introductions and removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14,
                 t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/R_r"+str(rg)+"_ir.txt",
               np.transpose(s1), fmt='%i')

############################################################################

# 2) 100 initial infected in five foci

S_T0 = np.array(XY['f5_t0'])

# a) No new events

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015)
    np.savetxt("./results/all_points/f5_r"+str(rg)+"_0.txt",
               np.transpose(s1), fmt='%i')

# b) Removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14)
    np.savetxt("./results/all_points/f5_r"+str(rg)+"_r.txt",
               np.transpose(s1), fmt='%i')

# c) New introductions

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/f5_r"+str(rg)+"_i.txt",
               np.transpose(s1), fmt='%i')

# d) New introductions and removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14,
                 t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/f5_r"+str(rg)+"_ir.txt",
               np.transpose(s1), fmt='%i')

############################################################################

# 3) 100 initial infected in one focus

S_T0 = np.array(XY['f1_t0'])

# a) No new events

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015)
    print(rg)
    np.savetxt("./results/all_points/f1_r"+str(rg)+"_0.txt",
               np.transpose(s1), fmt='%i')

# b) Removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/f1_r"+str(rg)+"_i.txt",
               np.transpose(s1), fmt='%i')

# c) New introductions

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14)
    np.savetxt("./results/all_points/f1_r"+str(rg)+"_r.txt",
               np.transpose(s1), fmt='%i')

# d) New introductions and removing individuals

for rg in RG_LIST:
    s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*20, t_save=6, beta=0.015,
                 rg=rg, nu=1, sigma=8, lmbda=0.015, t_remove=12*14,
                 t_new=12, n_new=N_NEW)
    np.savetxt("./results/all_points/f1_r"+str(rg)+"_ir.txt",
               np.transpose(s1), fmt='%i')

############################################################################
    