#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:06:22 2020

@author: suraj
"""

import numpy as np
np.random.seed(22)
import matplotlib.pyplot as plt

font = {'family' : 'Times New Roman',
        'size'   : 12}    
plt.rc('font', **font)

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

#%%
from mpi4py import MPI
comm = MPI.COMM_WORLD # create communicator
my_rank = comm.Get_rank() # get rank of current processor
p = comm.Get_size() # get totall number of processors

if my_rank != 0:
    message = 'Hello from '+str(my_rank)
    comm.send(message, dest=0)
else:
    for procid in range(1,p):
        message = comm.recv(source=procid)
        print('Process 0 receives message from process ', procid,':', message)
        
#%%
        