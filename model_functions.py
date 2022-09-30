# -*- coding: utf-8 -*-
"""
@author: Martina Cendoya

Functions for the disease spread simulation algorithm
"""

from random import uniform
import numpy as np
import numba as nb
import scipy.special as sc


def nb_dist(d_max, coords, indI, indS):
    """Neighbors within a distance radius of d_max"""
    nb_0 = []
    for i in indI:
        nb_i = []
        coords_i = coords[i]
        for j in indS:
            s = 0.0
            for k in range(2):
                tmp = coords_i[k] - coords[j, k]
                s += tmp * tmp
                dist = s**0.5
            if dist <= d_max:
                nb_i.append(j)
        nb_0.append(nb_i)
    return nb_0

@nb.njit(fastmath=True)
def dist_eu(coords_i, coords_j):
    """Euclidean distance"""
    s = 0.0
    for i in range(2):
        tmp = coords_i[i] - coords_j[i]
        s += tmp * tmp
    d = s**0.5
    return float(d)


@nb.njit(fastmath=True)
def cor_mat(d, rg, nu):
    """Matérn correlation function"""
    nu = float(nu)
    phi = (np.sqrt(8*nu)/rg)
    a = (2**(nu - 1) * sc.gamma(nu))**(-1)
    b = (d * phi)**nu * sc.kv(nu, d * phi)
    cor = a * b
    return cor

@nb.njit
def new_intro(st, nI):
    """New introductions"""
    S = 0
    I_new = 3
    if len(st[st == S]) > nI:
        idx = np.random.choice(np.where(st == S)[0], size=nI, replace=False)
    st[idx] = I_new
    return st

@nb.njit(parallel=True, fastmath=True)
def stat(s_init, coords, beta, rg, nu, sigma, lmbda, t_max, t_save, t_remove=None,
         t_new=None, n_new=None):
    """Function for spread simulation"""
    t_s = np.arange(0, t_max+1, t_save)
    t_s = t_s.astype(nb.int32)
    s_final = np.empty((len(t_s), len(s_init)), dtype=np.float64)
    s_final[0] = s_init
    s_t1 = np.copy(s_init)
    dist_max = rg + (rg/2)
    sint_expr = np.array([1]*len(coords))
    if t_remove is not None:
        t_inf = np.array([1]*len(coords))
    for t in range(1, t_max+1):
        s_t1[s_t1 == 3] = 1
        indI = np.where(s_t1 > 0)[0]
        indS = np.where(s_t1 == 0)[0]
        indIA = np.where(s_t1 == 1)[0]
        stat_nb = s_t1[s_t1 > 0]
        if t_remove is not None:
            for infected in nb.prange(len(indI)):
                if t_inf[indI[infected]] == t_remove:
                    s_t1[indI[infected]] = -1
                else:
                    t_inf[indI[infected]] += 1
        for i in nb.prange(len(indS)):
            dx = 0
            for j in range(len(indI)):
                dist = dist_eu(coords[indS[i]], coords[indI[j]])
                if dist < dist_max:
                    cor_m = cor_mat(dist, rg, nu)
                    if stat_nb[j] == 1:
                        dx += (cor_m * beta * lmbda)
                    else: dx += (cor_m * beta)
            dx_p = 1 - np.exp(-(dx))
            u = uniform(0, 1)
            if dx_p > u:
                s_t1[indS[i]] = 1
        for ia in nb.prange(len(indIA)):
            if sint_expr[indIA[ia]] == sigma:
                s_t1[indIA[ia]] = 2
            else:
                sint_expr[indIA[ia]] += 1
        if t_new is not None and t % t_new == 0:
            n_idx = int(t/t_new)-1
            nI = n_new[n_idx]
            s_t1 = new_intro(s_t1, nI)
        if t in t_s:
            idx = int(np.where(t_s == t)[0][0])
            print(idx, t)
            s_final[idx] = s_t1
    return s_final
