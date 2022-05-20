# MAIN SIMULATION LOOP

import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from input import *
from boundary import *
from cutoff import *
from files import *
from LJ import *
from neighbours import *
from parameters import *
from d_max import *

df = read_out("Initial_conf.out")

co = cutoff()

d_max = d_max()

for tt in range(N_step + 1):
    
    dff = df
    for i in range(N_atom):
        for j in range(3):
            dx  = (random.random() - 0.5) * d_max
            df.iloc[i,j] += dx
    
    if tt%10 == 0:
        sk = Skin(df,N_atom,co)
    
    neig   = Neighbours(df,sk,N_atom,co)
    
    if tt%10 == 0:
        skk  = Skin(dff,N_atom,co)
    
    neigh    = Neighbours(dff,skk,N_atom,co)
    
    U_x = lennard_jones(df,neig,N_atom)
    U   = lennard_jones(dff,neigh,N_atom)
    
    for i in range(N_atom):
        dU = U_x[i][0] - U[i][0]
        for j in range(3):
            if dU >= 0:
                p_acc = dU/kT
                if p_acc < math.sqrt(math.log(random.random())**2):
                    df.iloc[i,j] = dff.iloc[i,j]
    
    n_df = df
    periodic_cond(n_df,box,N_atom,bound)
    df = n_df
    
    if tt%N_write == 0:
        crds  = df.to_numpy()
        name  = str(tt) + '.out'
        np.savetxt(name, crds)   
        
#Author      : Kalith M Ismail.
#Objective   : Metropolis monte-carlo simulation of noble gas atoms in the box. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 20/05/2022.
