# -*- coding: utf-8 -*-
"""
Functions for spatial individual-based disease spread simulation.

Values for the disease states of individuals:
    - Susceptible individuals: S = 0.
    - Asymptomatic infected individuals: Ia = 0.015 (transmission rate reduction value).
    - Symptomatic infected individuals: Is = 1.
"""

import numpy as np
import numba as nb
import scipy.special as sc



@nb.njit(fastmath=True)
def dist_eu(coords_i, coords_j):
    """
    Euclidean distance function between two points ij.
    
    Args:
        - coords_i: numpy array of shape (1,2) with coordinates of point i.
        - coords_j: numpy array of shape (1,2) with coordinates of point j.
        
    Returns:
        - d: float Euclidean distance between two points.
    """
    s = 0.0
    for i in range(2):
        tmp = coords_i[i] - coords_j[i]
        s += tmp * tmp
    d = s**0.5
    return float(d)


def nb_dist(D_MAX, coords, ind, nbs):
    """
    Finds neighbors of individuals within a maximum distance of D_MAX.
    
    Args:
        - D_MAX: the maximum distance for two individuals to be considered neighbors.
        - coords: array of shape (n,2) containing the x,y coordinates of each individual.
        - ind: the index or indices of the individual(s) for which to find neighbors.
        - nbs: the index of all potential neighbors for the given individual(s).
        
    Returns:
        - nb_0: a list of lists, where the i-th list contains the indices of all neighbors for the i-th individual.
    """
    
    nb_0 = []
    for i in ind:
        nb_i = []
        for j in nbs:
            dist = dist_eu(coords[i], coords[j])
            if dist <= D_MAX:
                nb_i.append(j)
        nb_0.append(nb_i)
    return nb_0


@nb.njit(fastmath=True)
def cor_matern(dist, rg, nu):
    """
    Matérn correlation function between two points.
    
    Args:
        - dist: Euclidean distance between two points.
        - rg: range parameter.
        - nu: smoothing parameter.
        
    Returns:
        - cor: float spatial correlation between two points.
    """
    nu = float(nu)
    phi = (np.sqrt(8*nu)/rg)
    a = (2**(nu - 1) * sc.gamma(nu))**(-1)
    b = (dist * phi)**nu * sc.kv(nu, dist * phi)
    cor = a * b
    return cor


@nb.njit(fastmath=True)
def f_inf(coords, indS_i, indI, s_t, beta, rg, nu, dist_max):
    """
    Calculates the probability of infection for a susceptible individual based on the
    force of infection generated by infected individuals at distance < dist_max of the susceptible individual.
    
    Args:
        - coords: numpy array of shape (n,2) containing the x,y coordinates of each individual.
        - indS_i: integer index of the susceptible individual for which to calculate the probability of infection.
        - indI: integer array of indices of infected individuals.
        - s_t: numpy array of length n containing the status of each individual (0=S, 0.015=Ia, 1=Is).
        - beta: transmission rate parameter.
        - rg: range parameter of the Matérn correlation function.
        - nu: smoothing parameter of the Matérn correlation function.
        - dist_max: maximum distance to evaluate the correlation.
    
    Returns:
        - dx_p: float probability of infection for the susceptible individual.
    """
    
    # Initialize the force of infection to zero
    dx=0
    
    # Loop through all infected individuals
    for j in indI:
        # Calculate the Euclidean distance between the susceptible individual and the infected individual
        dist = dist_eu(coords[indS_i], coords[j])
        
        # If the distance is within the specified maximum distance
        if dist < dist_max:
            # Calculate the value of the Matérn correlation function
            cor_m = cor_matern(dist, rg, nu)
            # Adjust the transmission rate based on the status of the infected individual
            lmbda = s_t[j]
            # Calculate the force of infection generated by the infected individual
            dx += (cor_m * beta * lmbda)
            
    # Calculate the probability of infection based on the total force of infection
    dx_p = 1 - np.exp(-(dx))
    
    return dx_p

#
# Individual-based spread: probability of infection calculated for each susceptible individual
#

@nb.njit(parallel=True, fastmath=True)
def stat_t(s_t, coords, dist_max, rg, nu, beta, sigma, sint_expr):
    """
    Performs the spread of infection for a given time-step using the 'individual-based spread' model.
    
    Args:
        - s_t: numpy array with the status of each individual (0=S, 0.015=Ia, 1=Is).
        - coords: numpy array of shape (n, 2) with the x and y coordinates of each point.
        - dist_max: maximum distance to evaluate the correlation between individuals.
        - rg: range parameter of the Matérn correlation function.
        - nu: smoothing parameter of the Matérn correlation function.
        - beta: transmission rate parameter.
        - sigma: time from infection to symptom expression.
        - sint_expr: numpy array with the time since infection for each Ia individual.
        
    Returns:
        Two numpy arrays:
        - s_t: updated numpy array with the status of each individual after the spread.
        - sint_expr: updated numpy array with the time since infection for each Ia individual after the spread.
    """
    
    # Index of infected individuals (Ia and Is)
    indI = np.where(s_t > 0)[0]
    # Index of susceptible individuals (S)
    indS = np.where(s_t == 0)[0]
    # Index of asymptomatic infected individuals (Ia)
    indIA = np.where(s_t == 0.015)[0]
    
    # Loop over S individuals
    for i in nb.prange(len(indS)):
        # Probability of infection
        prob = f_inf(coords, indS[i], indI, s_t, beta, rg, nu, dist_max)
        # Random variable from a Bernoulli distribution with probability 'prob'
        u = np.random.binomial(1, prob)
        # If the Bernoulli random variable (u) is equal to 1, the individual becomes Ia
        s_t[indS[i]] = 0.015 if u==1 else 0
        
    # Loop over Ia individuals
    for j in nb.prange(len(indIA)):
        if sint_expr[indIA[j]] < sigma:
            # Update sint_expr if the individual has not yet expressed symptoms
            sint_expr[indIA[j]] += 1
        else:
            # The individual becomes symptomatic if the symptom expression time has passed
            s_t[indIA[j]] = 1
            
    # Return the updated status array and time since infection array
    return(s_t, sint_expr)


@nb.njit(parallel=True, fastmath=True)
def stat(s_init, coords, beta, rg, nu, sigma, t_max, t_save):
    """
    Run the 'individual-based spread' model from time t=0 to time t_max.
    
    Args:
        - s_init: numpy array of the initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
        - coords: numpy array of shape (n, 2) with the x and y coordinates of each point.
        - beta: transmission rate.
        - rg: range of the Matérn correlation function.
        - nu: smoothing parameter of the Matérn correlation function.
        - sigma: time from infection to symptom expression.
        - t_max: maximum time for spread.
        - t_save: interval between time-steps to save the status of individuals.
        
    Returns:
        - s_final: array containing the status of individuals at different time-steps.
    """
    
    # Create array of time-steps to save
    t_s = np.arange(0, t_max+1, t_save)
    t_s = t_s.astype(nb.int32)
    
    # Initialize array to store final results
    s_final = np.empty((len(t_s), len(s_init)), dtype=np.float64)
    s_final[0] = s_init
    
    # Copy initial status array for processing
    s_t = np.copy(s_init)
    
    # Set maximum distance to evaluate correlation
    dist_max = rg + (rg/2)
    
    # Initialize array to store time since infection
    sint_expr = np.array([0]*len(coords))
    sint_expr[s_init>0] = 1
    
    # Loop over time-steps
    for t in range(1, t_max+1):
        # Run the 'stat_t' function to spread the disease for one time-step
        s_t, sint_expr = stat_t(s_t, coords, dist_max, rg, nu, beta, sigma, sint_expr)
        
        # Save the status array if current time-step is in the list of time-steps to save
        if t in t_s:
            idx = int(np.where(t_s == t)[0][0])
            s_final[idx] = s_t
            
    # Return array containing the status of individuals at different time-steps
    return s_final

#
# Cell-based model: probability of infection calculated for each cell
#

@nb.njit(parallel=True, fastmath=True)
def stat_tc(s_t, idx_cells, coords, dist_max, rg, nu, beta, sigma, sint_expr):
    """
    Performs the spread of infection for a given time-step using the optimized 'cell-based spread' model.
    
    Args:
        - s_t: numpy array with the status of each individual (0=S, 0.015=Ia, 1=Is).
        - idx_cells: numpy array with the cell index of each individual.
        - coords: numpy array of shape (n, 2) with the x and y coordinates of each point.
        - dist_max: maximum distance to evaluate the correlation between individuals.
        - rg: range parameter of the Matérn correlation function.
        - nu: smoothing parameter of the Matérn correlation function.
        - beta: transmission rate parameter.
        - sigma: time from infection to symptom expression.
        - sint_expr: numpy array with the time since infection for each Ia individual.
    
    Returns:
        Two numpy arrays:
        - s_t: updated numpy array with the status of each individual after the spread.
        - sint_expr: updated numpy array with the time since infection for each Ia individual after the spread.
    """
    
    # Index of infected individuals (Ia and Is)
    indI = np.where(s_t > 0)[0]
    # Index of susceptible individuals (S)
    indS = np.where(s_t == 0)[0]
    # Index of asymptomatic infected individuals (Ia)
    indIA = np.where(s_t == 0.015)[0]
    # cell index of each S individual
    indS_cell = idx_cells[indS]
    # Probability of infection per cell
    cells = np.unique(indS_cell)
    
    # Loop over cells
    for cell in nb.prange(len(cells)):
        # Random selection of an S individual in the cell
        i_cell = np.random.choice(indS[indS_cell == cells[cell]], 1)[0]
        # Probability of infection for the selected individual
        prob = f_inf(coords, i_cell, indI, s_t, beta, rg, nu, dist_max)
        
        # Loop over S individuals in the cell
        for i in indS[indS_cell == cells[cell]]:
            # Random variable from a Bernoulli distribution with probability 'prob'
            u = np.random.binomial(1, prob)
            # If the Bernoulli random variable (u) is equal to 1, the individual becomes Ia
            s_t[i] = 0.015 if u==1 else 0
            
    # Loop over Ia individuals
    for j in nb.prange(len(indIA)):
        if sint_expr[indIA[j]] < sigma:
            # Update sint_expr if the individual has not yet expressed symptoms
            sint_expr[indIA[j]] += 1
        else:
            # The individual becomes symptomatic if the symptom expression time has passed
            s_t[indIA[j]] = 1
            
    # Return the updated status array and time since infection array
    return(s_t, sint_expr)


@nb.njit(parallel=True, fastmath=True)
def stat_c(s_init, idx_cells, coords, beta, rg, nu, sigma, t_max, t_save):
    """
    Run the 'cell-based spread' model from time t=0 to time t_max.
    
    Args:
        - s_init: array of the initial status of individuals at t=0 (0 = S, 0.015 = Ia, 1 = Is).
        - idx_cells: array with the cell index of each individual.
        - coords: array of coordinates of points.
        - beta: transmission rate.
        - rg: range of the Matérn correlation function.
        - nu: smoothing parameter of the Matérn correlation function.
        - sigma: time from infection to symptom expression.
        - t_max: maximum time for spread.
        - t_save: interval between time-steps to save the status of individuals.
        
    Returns:
        - s_final: array containing the status of individuals at different time-steps.
    """
    # Create array of time-steps to save
    t_s = np.arange(0, t_max+1, t_save)
    t_s = t_s.astype(nb.int32)
    
    # Initialize array to store final results
    s_final = np.empty((len(t_s), len(s_init)), dtype=np.float64)
    s_final[0] = s_init
    
    # Copy initial status array for processing
    s_t = np.copy(s_init)
    
    # Set maximum distance to evaluate correlation
    dist_max = rg + (rg/2)
    
    # Initialize array to store time since infection
    sint_expr = np.array([0]*len(coords))
    sint_expr[s_init>0] = 1
    
    # Loop over time-steps
    for t in range(1, t_max+1):
        # Run the 'stat_tc' function to spread the disease for one time-step
        s_t, sint_expr = stat_tc(s_t, idx_cells, coords, dist_max, rg, nu, beta, sigma, sint_expr)
        
        # Save the status array if current time-step is in the list of time-steps to save
        if t in t_s:
            idx = int(np.where(t_s == t)[0][0])
            s_final[idx] = s_t
            
    # Return array containing the status of individuals at different time-steps
    return s_final
