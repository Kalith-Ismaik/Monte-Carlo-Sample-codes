#Prediction of d_max value for metropolis monte carlo simulation.

import math
import random
import matplotlib.pyplot as plt

from LJ import *

def d_max():
	acc   = []
	dist  = []
	d     = 2.8329999999999997
	d_max = 1
	U     = U_pot(d,lj_s,lj_e)

	while (d_max < 101):
    		per = 0
    		rand = random.random()
    		for i in range(10000):
        		dx  = (random.random() - 0.5) * d_max
        		dx += d
		        U_x = U_pot(dx,lj_s,lj_e)
		        dU  = U_x - U
		        if dU < 0:
		            per += 1
		        elif dU >= 0:
		            p_acc = dU/(kT)
		            if p_acc > math.sqrt(math.log(random.random())**2):
	 	               per += 1
    		per /= 100
    		#print("Displaced_distance:",d_max,"  Acceptance_rate:",per,"%")
    		acc.append(per)
    		dist.append(d_max)
    		d_max += 1

	loc     = acc.index(max(acc))
	d_max   = dist[loc]
	print("d_max value the sample crystal is considered to be:",dist[loc],"A")
	print("Acceptance rate of the considered d_max:",acc[loc],"%")
	print("d_max is predicted with the nearest value to the cutoff:","50","%")
	fig, ax = plt.subplots()
	ax.plot(dist,acc)
	plt.xlabel("Displaced_distance [A]")
	plt.ylabel("Acceptance_rate [%]")
	plt.title("Acceptance rate vs Displaced distance")
	plt.show()
	
	return d_max

#Author      : Kalith M Ismail.
#Objective   : Prediction of metropolis monte carlo simulation.
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 20/05/2022.
