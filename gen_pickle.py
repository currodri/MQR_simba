#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 11:14:07 2019

@author: currorodriguez
"""
# Import required libraries
import numpy as np
import pickle
import sys

MODEL = sys.argv[1]  # e.g. m50n512

progen_file = '../progen_analysis/%s/progen_%s.pkl' % (MODEL, MODEL)

# Extract progen data from txt files
obj = open(progen_file, 'rb')
d = pickle.load(obj)
ngal = d['galaxies_per_snap'][0]
print('Total number of galaxies at z = 0: '+str(ngal))

#Store the galaxies sorted in objects of type GalaxyData
galaxies = []
for i in range(ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['t_gal' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]
    gal_type = d['gal_type'+str(i)][::-1]
    gal_pos = d['gal_pos'+str(i)][::-1]
    caesar_id = d['caesar_id'+str(i)][::-1]
    galaxy = GalaxyData(i, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type, gal_pos, caesar_id)
    galaxies.append(galaxy)

max_ngal = len(galaxies)
mass_limit = 9.5
min_merger_ratio = 0.2
max_redshift_mergers = 2.5

# Perform the search for mergers
mergers, sf_galaxies = merger_finder(galaxies, min_merger_ratio, 10**mass_limit, max_redshift_mergers, out_file=True)

print('Merger analysis done.')

# Perform the quenching and rejuvenation analysis
galaxies_interpolated = quenchingFinder2(galaxies[0:max_ngal], 1, mass_limit)

quenchingFinder2(galaxies_interpolated, 1, mass_limit, interpolation=True, out_file=True)

print('Quenching analysis done.')
