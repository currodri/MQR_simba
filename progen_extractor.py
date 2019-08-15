#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

This code is the basis of the evolution analysis of major quantities in the SIMBA simulation. Given the relation of galaxies
at z = 0 with their progenitors in previous snapshots, the information coming from the Caesar files can be linked together.

All this information is saved in a dictionary that is dumped into a pickle file.

@author: currorodriguez
"""

# Import required libraries
import caesar
import numpy as np
from functools import reduce
import os
from astropy.cosmology import FlatLambdaCDM
import cPickle as pickle
import sys

MODEL = sys.argv[1]  # e.g. m50n512
WIND = sys.argv[2]  # e.g. s50 for Simba

caesarfile = '/home/rad/data/%s/%s/Groups/' % (MODEL,WIND)
progenref_file = '/disk01/rad/sim/%s/%s/Groups/progen_%s_151.dat' % (MODEL,WIND,MODEL)
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../progen_analysis/%s/' % (MODEL)


progenref = open(progenref_file, 'r').readlines()
lengal = int(progenref[0].split(' ')[0])
progenref_data = []
lines = int((len(progenref)-1)/(2*lengal))
print(lengal)
for galaxy in range(0, lengal):

    start = galaxy*lines + 1
    end = start + lines

    super_line = reduce(lambda x,y:x+y,[progenref[line] for line in range(start, end)])
    super_line = super_line.replace('[', ' ')
    super_line = super_line.replace(']', ' ')

    super_line = super_line.split()

    progenref_data.append([int(x) for x in super_line])

d = {}
for j in range(0, lengal):
    d['m'+str(j)] = np.array([])
    d['sfr'+str(j)] = np.array([])
    d['h1_gas'+str(j)] = np.array([])
    d['h2_gas'+str(j)] = np.array([])
    d['g_type'+str(j)] = np.array([])
    d['caesar_id'+str(j)] = np.array([])
    d['bhm'+str(j)] = np.array([])
    d['bhar'+str(j)] = np.array([])
    d['local_den'+str(j)] = np.array([])
    d['z'+str(j)] = np.array([])
    d['t'+str(j)] = np.array([])


snaps = filter(lambda file:file[-5:]=='.hdf5' and file[0]=='m' and int(file[-8:-5])!=116, os.listdir(caesarfile))
snaps_sorted = sorted(snaps,key=lambda file: int(file[-8:-5]), reverse=True)
print('Progenitor indexes obtained from .dat file.')
print('Saving data to dictionary...')
d['sf_galaxies_per_snap'] = np.zeros(len(snaps_sorted))
d['galaxies_per_snap'] = np.zeros(len(snaps_sorted))
d['sf_galaxies_mass'] = np.array([])
d['redshifts'] = np.zeros(len(snaps_sorted))
d['t_hubble'] = np.zeros(len(snaps_sorted))
print(snaps_sorted)
def sfr_condition(type, time):
    if type == 'start':
        lsfr = np.log10(1/(time))-9
    elif type == 'end':
        lsfr  = np.log10(0.2/(time))-9
    return lsfr

for s in range(0, len(progenref_data[0])+1):

    #if snaps_sorted[s] != 'm50n512_116.hdf5' and WIND == 's50': # This condition is just to solve an issue with that snapshot
    sim = caesar.load(caesarfile+snaps_sorted[s],LoadHalo=False) # load caesar file

    # initialize simulation parameters
    redshift = sim.simulation.redshift  # this is the redshift of the simulation output
    h = sim.simulation.hubble_constant  # this is the hubble parameter = H0/100
    cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)  # set our cosmological parameters
    thubble = cosmo.age(redshift).value  # age of universe at this redshift

    # Get galaxy info
    gals = np.asarray([i.central for i in sim.galaxies])   # read in galaxies from caesar file
    ms = np.asarray([i.masses['stellar'] for i in sim.galaxies])   # read in stellar masses of galaxies
    mHI = np.asarray([i.masses['HI'] for i in sim.galaxies])   # read in neutral hydrogen masses
    mH2 = np.asarray([i.masses['H2'] for i in sim.galaxies])   # read in molecular hydrogen
    sfr = np.asarray([i.sfr for i in sim.galaxies])   # read in instantaneous star formation rates
    galpos = np.array([g.pos.d for g in sim.galaxies]) # the .d removes the units
    caesar_id = np.array([i.GroupID for i in sim.galaxies]) # getting the Caesar ID for each galaxy
    local_mass_density = np.array([i.local_mass_density for i in sim.galaxies]) # getting environmental measure of mass density
    bh_dot = [] # getting the BH accretion rate for the most massive BH particle
    bh_mass = [] # getting mass of most massive BH particle in galaxy
    for gal in sim.galaxies:
        try:
            bh_dot.append(float(gal.bhmdot.d))
        except AttributeError:
            bh_dot.append(0.0)
    for gal in sim.galaxies:
        try:
            bh_mass.append(float(gal.masses['bh'].d))
        except KeyError:
            bh_mass.append(0.0)
    bh_dot = np.asarray(bh_dot)
    bh_mass = np.asarray(bh_mass)
    ssfr_gal = sfr/ms
    sfgals = 0
    ssfr_cond = sfr_condition('end',thubble)
    sfmass = []
    for i in range(0, len(gals)):
        if ssfr_gal[i] >= 10**ssfr_cond:
            sfgals = sfgals + 1
            sfmass.append(ms[i])
    sfmass = np.asarray(sfmass)
    d['sf_galaxies_mass'] = np.concatenate((d['sf_galaxies_mass'],sfmass))
    print('Number of star forming galaxies in this snapshot: '+str(sfgals))
    print('Median mass of star forming galaxies in this snapshot: '+str(np.median(d['sf_galaxies_mass'][s]))+' M*')
    d['sf_galaxies_per_snap'][s] = sfgals
    d['galaxies_per_snap'][s] = len(gals)
    d['redshifts'][s] = redshift
    d['t_hubble'][s] = thubble
    if s==0:
        d['boxsize_in_kpccm'] = sim.simulation.boxsize.to('kpccm')
    for k in range(0, lengal):
        if s==0:
            d['m'+str(k)] = np.concatenate((d['m'+str(k)], ms[k]), axis=None)
            d['sfr'+str(k)] = np.concatenate((d['sfr'+str(k)],sfr[k]), axis=None)
            d['h1_gas'+str(k)] = np.concatenate((d['h1_gas'+str(k)],mHI[k]), axis=None)
            d['h2_gas'+str(k)] = np.concatenate((d['h2_gas'+str(k)],mH2[k]), axis=None)
            d['z'+str(k)] = np.concatenate((d['z'+str(k)],redshift), axis=None)
            d['t'+str(k)] = np.concatenate((d['t'+str(k)],thubble), axis=None)
            d['g_type'+str(k)] = np.concatenate((d['g_type'+str(k)],gals[k]), axis=None)
            d['caesar_id'+str(k)] = np.concatenate((d['caesar_id'+str(k)],caesar_id[k]), axis=None)
            d['bhar'+str(k)] = np.concatenate((d['bhar'+str(k)],bh_dot[k]), axis=None)
            d['bhm'+str(k)] = np.concatenate((d['bhm'+str(k)],bh_mass[k]), axis=None)
            d['pos'+str(k)] = np.array([galpos[k]])
            d['local_den'+str(k)] = np.concatenate((d['local_den'+str(k)],local_mass_density[k]), axis=None)

        elif progenref_data[k][s-1] !=-1:
            index = progenref_data[k][s-1]
            d['m'+str(k)] = np.concatenate((d['m'+str(k)], ms[index]), axis=None)
            d['sfr'+str(k)] = np.concatenate((d['sfr'+str(k)],sfr[index]), axis=None)
            d['h1_gas'+str(k)] = np.concatenate((d['h1_gas'+str(k)],mHI[index]), axis=None)
            d['h2_gas'+str(k)] = np.concatenate((d['h2_gas'+str(k)],mH2[index]), axis=None)
            d['z'+str(k)] = np.concatenate((d['z'+str(k)],redshift), axis=None)
            d['t'+str(k)] = np.concatenate((d['t'+str(k)],thubble), axis=None)
            d['g_type'+str(k)] = np.concatenate((d['g_type'+str(k)],gals[index]), axis=None)
            d['caesar_id'+str(k)] = np.concatenate((d['caesar_id'+str(k)],caesar_id[index]), axis=None)
            d['bhar'+str(k)] = np.concatenate((d['bhar'+str(k)],bh_dot[index]), axis=None)
            d['bhm'+str(k)] = np.concatenate((d['bhm'+str(k)],bh_mass[index]), axis=None)
            d['pos'+str(k)] = np.concatenate((d['pos'+str(k)],np.asarray([galpos[index]])), axis=0)
            d['local_den'+str(k)] = np.concatenate((d['local_den'+str(k)],local_mass_density[index]), axis=None)
print('Data saved to dictionary.')
output = open(results_folder+'progen_'+str(MODEL)+'.pkl','wb')
pickle.dump(d, output)
print('Data saved in pickle file.')
output.close()
print('Progen extraction of galactic data: DONE!')
