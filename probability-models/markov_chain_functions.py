# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 08:03:10 2022

@author: sequo
"""
import numpy as np

def long_run_proportion(p):
    b = np.zeros(np.shape(p)[0])
    b[-1] = 1
    I = np.identity(np.shape(p)[0])
    a = p.T - I
    a[-1] = np.array([1 for i in range(np.shape(p)[0])])
    pi = np.linalg.inv(a).dot(b)
    return pi

def gamblers_ruin(p, i, n):
    q = 1-p
    num = 1 - ((q/p)**i)
    den = 1 - ((q/p)**n)
    return num/den

def e_x_mc(p, n, a):
    p_n = np.linalg.matrix_power(p, n)
    p_n_x = np.zeros(len(p_n))
    for i in range(len(p_n)):
        p_n_x[i] = p_n[:,i].dot(a)
    e_x = p_n_x.dot(np.array([i for i in range(len(p_n))]))
    return e_x

def calc_transient_times(pt):
    s = np.linalg.inv(np.identity(np.shape(pt)[0])-pt)
    return s

def calc_transient_total_time(pt):
    s = calc_transient_times(pt)
    g = s.dot(np.ones(np.shape(pt[0])))
    return g

def calc_exit_distribution(pt, a):
    s = calc_transient_times(pt)
    h = s.dot(a)
    return h

def check_time_reversible(p):
    pi = long_run_proportion(p)
    D = np.diag(pi)
    left = D.dot(p).round(3)
    right = p.T.dot(D).round(3)
    return np.array_equal(left, right), left, right
