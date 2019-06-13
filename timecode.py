!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 13 Jun 2019

@author: currorodriguez
"""

# Import required libraries
import caesar
import numpy as np
from functools import reduce
import os
from astropy.cosmology import FlatLambdaCDM

caesarfile = '/home/rad/data/m100n1024/s50/Groups/' #input('Final group file: ')
snaps = filter(lambda file:file[-5:]=='.hdf5' and file[0]=='m', os.listdir(caesarfile))
snaps_sorted = sorted(snaps,key=lambda file: int(file[-8:-5]), reverse=True)

timefile = open('times_m100n1024.txt','w')
redfile = open('redshifts_m100n1024.txt','w')

for s in range(0, len(snaps_sorted)):
    sim = caesar.load(caesarfile+snaps_sorted[s],LoadHalo=False) # load caesar file
    # initialize simulation parameters
    redshift = sim.simulation.redshift  # this is the redshift of the simulation output
    cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)  # set our cosmological parameters
    thubble = cosmo.age(redshift).value
    timefile.write(str(thubble)+'\n')
    redfile.write(str(redshift)+'\n')
timefile.close()
redfile.close()
