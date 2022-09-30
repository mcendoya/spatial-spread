# -*- coding: utf-8 -*-
"""
@author: Martina Cendoya

Simulations in the test area
"""

from random import sample
import numpy as np
import pandas as pd
from numba import set_num_threads
import model_functions as fn

set_num_threads(4)

#######################################################################
##########################   TEST AREA   ##############################
#######################################################################

XY = pd.read_csv("./Data/coords_ex.csv")
COORDS = XY[['x', 'y']].values

RG_LIST = [100, 200, 500]
NU_LIST = [0.5, 1, 1.5]
BETA_LIST = [0.008, 0.015, 0.04]

########################################################################
# 1) PARAMETERS
########################################################################

# 1.1) 10 random initial infected

S_T0 = np.array(XY['R_t0'])

# a) Monthly simulations

for beta in BETA_LIST:
    for rg in RG_LIST:
        for nu in NU_LIST:
            s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*10, t_save=6,
                         beta=beta, rg=rg, nu=nu, sigma=8, lmbda=0.015)
            np.savetxt("./results/area_ex/parameters/R_b"+str(beta)+
                       "_r"+str(rg)+"_nu"+str(nu)+".txt",
                       np.transpose(s1), fmt='%i')

# b) Daily simulations

for rg in RG_LIST:
    for nu in NU_LIST:
        s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=365*10, t_save=182.5,
                     beta=0.00036, rg=rg, nu=nu, sigma=242, lmbda=0.015)
        np.savetxt("./results/area_ex/parameters/R_b0.00036_r"+
                   str(rg)+"_nu"+str(nu)+".txt",
                   np.transpose(s1), fmt='%i')

###########################################################################

# 1.2) 10 initial infected in one focus

S_T0 = np.array(XY['f1_t0'])

# a) Monthly simulations

for beta in BETA_LIST:
    for rg in RG_LIST:
        for nu in NU_LIST:
            s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=12*10, t_save=6,
                         beta=beta, rg=rg, nu=nu, sigma=8, lmbda=0.015)
            np.savetxt("./results/area_ex/parameters/f1_b"+str(beta)+
                       "_r"+str(rg)+"_nu"+str(nu)+".txt",
                       np.transpose(s1), fmt='%i')

# b) Daily simulations

for rg in RG_LIST:
    for nu in NU_LIST:
        s1 = fn.stat(s_init=S_T0, coords=COORDS, t_max=365*10, t_save=182.5,
                     beta=0.00036, rg=rg, nu=nu, sigma=242, lmbda=0.015)
        np.savetxt("./results/area_ex/parameters/f1_b0.00036_r"+
                   str(rg)+"_nu"+str(nu)+".txt",
                   np.transpose(s1), fmt='%i')

########################################################################
# 2) VARIABILITY
########################################################################

# 2.1) 10 random initial infected

for rg in RG_LIST:
    it = 0
    while it < 100:
        print("iteration = ", it+1)
        s_t0 = np.zeros(len(XY))
        idx = sample(range(0, len(s_t0)), 10)
        s_t0[idx] = 1
        s1 = fn.stat(s_init=s_t0, coords=COORDS, t_max=12*10, t_save=6,
                     beta=0.015, rg=rg, nu=1, sigma=8, lmbda=0.015)
        np.savetxt("./results/area_ex/variability/R_r"+str(rg)+
                   "s"+str(it+1)+".txt",
                   np.transpose(s1), fmt='%i')
        it += 1

##########################################################################

# 2.2) 10 initial infected in one focus

D_MAX = 1000 # Focus radious (m)

for rg in RG_LIST:
    it = 0
    while it < 100:
        print("iteration = ", it+1)
        s_t0 = np.zeros(len(XY))
        idx = sample(range(0, len(s_t0)), 1)
        s_t0[idx] = 1
        indI = idx
        indS = np.where(s_t0 == 0)[0]
        nb_1000 = fn.nb_dist(D_MAX, COORDS, indI, indS)
        nb_I = sample(nb_1000[0], 9)
        s_t0[nb_I] = 1
        s1 = fn.stat(s_init=s_t0, coords=COORDS, t_max=12*10, t_save=6,
                     beta=0.015, rg=rg, nu=1, sigma=8, lmbda=0.015)
        np.savetxt("./results/area_ex/variability/f1_r"+
                   str(rg)+"s"+str(it+1)+".txt",
                   np.transpose(s1), fmt='%i')
        it += 1

##########################################################################
