# -*- coding: utf-8 -*-
"""
Perform simulations of the spatial individual-based model

Parameters:
    - beta: transmission rate
    - rg: range of the Matérn correlation function.
    - nu: smoothing parameter of the Matérn correlation function.
    - sigma: time from infection to symptom expression.
"""

import functions_spread as fn
import pandas as pd
import numpy as np
from numba import set_num_threads

# Set number of threads to run in parallel
set_num_threads(4) 

# # Data with coordinates of individuals, cell index and initial introductions 
xy = pd.read_csv("coords.csv")
# Coordinates
coords = xy[['x', 'y']].values
# Cell index
idx_cells = np.array(xy['cell_ha'])


# t_max: maximum time for spread.
# t_save: interval between time-steps to save the status of individuals.
t_max, t_save = 30*12, 6

"""
PARAMETER EFFECTS

Perform a simulation with the optimized algorithm (probability of infection per cell), 
for different values of the parameters for each type of initial introduction (random, 5 foci and one focus).
"""

# Directory to store
dir_s = "./results/parameters/"

# Parameters
beta_val = [0.005, 0.015, 0.03]
rg_val = [250, 500, 750, 1000]
nu_val = [0.5, 1, 1.5]
sigma = 8

# RANDOM
ii = "R"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["R_1"])
s_init[s_init==1] = 0.015

# Run for each combination of the parameter values
for beta in beta_val:
    for rg in rg_val:
        for nu in nu_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_nu"+str(nu)+".txt", np.transpose(s1), fmt='%f')


# 5 FOCI
ii = "f5"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f5_1"])
s_init[s_init==1] = 0.015

# Run for each combination of the parameter values
for beta in beta_val:
    for rg in rg_val:
        for nu in nu_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_nu"+str(nu)+".txt", np.transpose(s1), fmt='%f')
            
# 1 FOCUS
ii = "f1"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f1_1"])
s_init[s_init==1] = 0.015

# Run for each combination of the parameter values
for beta in beta_val:
    for rg in rg_val:
        for nu in nu_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_nu"+str(nu)+".txt", np.transpose(s1), fmt='%f')

"""
INTRODUCTION VARIABILITY

Simulations with the optimized algorithm (probability of infection per cell), 
for different parameter values for each type of initial introduction (random, 5 foci and one focus). 
In each simulation different location of the initial infected individuals.
"""

# Directory to store
dir_s = "./results/introductions_10/"

# Parameters
beta_val = [0.005, 0.015, 0.03]
rg_val = [250, 500, 750, 1000]
nu = 1
sigma = 8

# RANDOM
ii = "R"

for x in range(1,11):
    # Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
    s_init = np.array(xy["R_{0}".format(x)])
    s_init[s_init==1] = 0.015
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_"+str(x)+".txt", np.transpose(s1), fmt='%f')
            
# 5 FOCI
ii = "f5"

for x in range(1,11):
    # Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
    s_init = np.array(xy["f5_{0}".format(x)])
    s_init[s_init==1] = 0.015
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_"+str(x)+".txt", np.transpose(s1), fmt='%f')
            
# 1 FOCUS
ii = "f1"

for x in range(1,11):
    # Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
    s_init = np.array(xy["f1_{0}".format(x)])
    s_init[s_init==1] = 0.015
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_"+str(x)+".txt", np.transpose(s1), fmt='%f')

"""
INTRINSIC VARIABILITY

Simulations with the optimized algorithm (probability of infection per cell), 
for different parameter values for each type of initial introduction (random, 5 foci and one focus). 
10 runs.
"""

# Directory to store
dir_s = "./results/simulations_10/"

# Parameters
beta_val = [0.005, 0.015, 0.03]
rg_val = [250, 500, 750, 1000]
nu = 1
sigma = 8

# RANDOM
ii = "R"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["R_1"])
s_init[s_init==1] = 0.015


for x in range(1,11):
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_s"+str(x)+".txt", np.transpose(s1), fmt='%f')

# 5 FOCI
ii = "f5"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f5_1"])
s_init[s_init==1] = 0.015


for x in range(1,11):
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_s"+str(x)+".txt", np.transpose(s1), fmt='%f')
            
# 1 FOCUS
ii = "f1"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f1_1"])
s_init[s_init==1] = 0.015


for x in range(1,11):
    for beta in beta_val:
        for rg in rg_val:
            # run function
            s1 = fn.stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save)
            # save results
            np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+"_s"+str(x)+".txt", np.transpose(s1), fmt='%f')


"""
INDIVIDUAL PROBABILITY (OPTIMIZATION VALIDATION)

Simulations with the algorithm that computes the probability of infection for 
each susceptible individual (non-optimized).
"""

# Directory to store
dir_s = "./results/individual/"

# Parameters
beta = 0.03
rg = 1000
nu = 1
sigma = 8

# RANDOM
ii = "R"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["R_1"])
s_init[s_init==1] = 0.015

# run function
s1 = fn.stat(s_init, coords, beta, rg, nu, sigma, t_max, t_save)
# save results
np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+".txt", np.transpose(s1), fmt='%f')


# 5 FOCI
ii = "f5"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f5_1"])
s_init[s_init==1] = 0.015

# run function
s1 = fn.stat(s_init, coords, beta, rg, nu, sigma, t_max, t_save)
# save results
np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+".txt", np.transpose(s1), fmt='%f')


# 1 FOCI
ii = "f1"
# Initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
s_init = np.array(xy["f1_1"])
s_init[s_init==1] = 0.015

# run function
s1 = fn.stat(s_init, coords, beta, rg, nu, sigma, t_max, t_save)
# save results
np.savetxt(dir_s+ii+"_r"+str(rg)+"_b"+str(beta)+".txt", np.transpose(s1), fmt='%f')
