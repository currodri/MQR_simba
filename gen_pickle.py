#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: Curro Rodriguez Montero, School of Physics and Astronomy,
            University of Edinburgh, JCMB, King's Buildings

This code uses the mergerFinder and quenchingFinder algorithms to create dictionaries of GalaxyData holding 
the details of mergers and quenching galaxies found. This dictionaries are saved in pickle file such that they
can be used afterwards mutiple times.

For questions about the code:
s1650043@ed.ac.uk
"""
# Import required libraries
import numpy as np
import pickle
import sys

from galaxy_class import GalaxyData
from mergerFinder import merger_finder
from quenchingFinder import quenchingFinder

MODEL = sys.argv[1]  # e.g. m50n512

progen_file = '../progen_analysis/%s/progen_%s.pkl' % (MODEL, MODEL) # File holding the progen info of galaxies

# Extract progen data from txt files
obj = open(progen_file, 'rb')
d = pickle.load(obj)
ngal = int(d['galaxies_per_snap'][0])
print('Total number of galaxies at z = 0: '+str(ngal))

#Store the galaxies sorted in objects of type GalaxyData
d_results = {}
d_results['redshifts'] = d['redshifts']
d_results['sf_galaxies_mass'] = d['sf_galaxies_mass']
d_results['sf_galaxies_per_snap'] = d['sf_galaxies_per_snap']
d_results['galaxies'] = []
for i in range(0,ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['t_gal' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]
    gal_type = d['gal_type'+str(i)][::-1]
    gal_pos = d['gal_pos'+str(i)][::-1]
    caesar_id = d['caesar_id'+str(i)][::-1]
    h1_gas = 0
    h2_gas = d['h2_gal'+str(i)][::-1] * d['m_gal'+str(i)][::-1]
    bh_m = 0
    bhar = 0
    galaxy = GalaxyData(i, sfr_gal, galaxy_m, z_gal, galaxy_t, h1_gas, h2_gas, bh_m, bhar, gal_type, gal_pos, caesar_id)
    d_results['galaxies'].append(galaxy)

# Setting the limiting conditions of the survey
max_ngal = len(d_results['galaxies'])
mass_limit = 9.5
min_merger_ratio = 0.2
max_redshift_mergers = 2.5

d_results['mass_limit'], d_results['min_merger_ratio'], d_results['max_redshift_mergers'] = mass_limit,min_merger_ratio,max_redshift_mergers

# Perform the search for mergers
d_results['galaxies'] = merger_finder(d_results['galaxies'][0:max_ngal], min_merger_ratio, 10**mass_limit, max_redshift_mergers)

print('Merger analysis done.')

# Perform the quenching and rejuvenation analysis
d_results['galaxies'] = quenchingFinder(d_results['galaxies'][0:max_ngal], 1, mass_limit)

print('Performing interpolation of quenching data...')

d_results['galaxies'] = quenchingFinder(d_results['galaxies'][0:max_ngal], 1, mass_limit, interpolation=True)

print('Quenching analysis done.')

print('Now, saving data to pickle file...')

import cPickle as pickle

output = open('./mandq_results_'+str(MODEL)+'.pkl','wb')
pickle.dump(d_results, output)
print('Data saved in pickle file.')
output.close()
