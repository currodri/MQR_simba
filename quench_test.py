#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 6 18:56:07 2019

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import pickle
import sys
import matplotlib
from scipy import stats
from astropy.cosmology import FlatLambdaCDM
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

from galaxy_class import GalaxyData
from mergerFinder import merger_finder
from quenchingFinder import quenchingFinder

MODEL = sys.argv[1]  # e.g. m50n512
WIND = sys.argv[2]   # e.g. s50
GALAXY = int(sys.argv[3])  # e.6. 4973

# Extract progen data from txt files
progen_file = '../progen_analysis/%s/progen_%s.pkl' % (MODEL, MODEL) # File holding the progen info of galaxies
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

sfr_gal = d['sfr' + str(GALAXY)][::-1]
z_gal = d['z' + str(GALAXY)][::-1]
galaxy_t = d['t' + str(GALAXY)][::-1]
galaxy_m = d['m'+str(GALAXY)][::-1]
gal_type = d['g_type'+str(GALAXY)][::-1]
gal_pos = d['pos'+str(GALAXY)][::-1]
caesar_id = d['caesar_id'+str(GALAXY)][::-1]
h1_gas = d['h1_gas'+str(GALAXY)][::-1]
h2_gas = d['h2_gas'+str(GALAXY)][::-1]
local_den = d['local_den'+str(GALAXY)][::-1]
bh_m = d['bhm'+str(GALAXY)][::-1]
bhar = d['bhar'+str(GALAXY)][::-1]
galaxy = GalaxyData(GALAXY, sfr_gal, galaxy_m, z_gal, galaxy_t, h1_gas, h2_gas, bh_m, bhar,local_den, gal_type, gal_pos, caesar_id)
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


galaxy = d_results['galaxies'][0]

above = np.zeros(len(galaxy.t[0]))
below = np.zeros(len(galaxy.t[0]))

for i in range(0, len(galaxy.t[0])):
	above[i] = np.log10(1/(galaxy.t[0][i]))-9
	below[i] = np.log10(0.2/(galaxy.t[0][i])) - 9


# Plot the results

fig, axes = plt.subplots(2, 1, sharex='col', figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
axes[0].plot(galaxy.t[0], np.log10(galaxy.ssfr[0]+1e-14), 'k-')
axes[0].plot(galaxy.t[0], above, 'b--', label=r'Star-forming threshold: sSFR $=1/t_{U}$')
axes[0].plot(galaxy.t[0], below, 'r--', label=r'Quench threshold: sSFR $=0.2/t_{U}$')
axes[0].plot(galaxy.t[1], np.log10(galaxy.ssfr[1]+1e-14), linestyle='--', color='grey', alpha=0.7)
mergers_idx = np.array([i.indx for i in galaxy.mergers])
reju_id = np.array([i for i in galaxy.rejuvenations])
thubble_start = np.array([galaxy.t[1][i.above9] for i in galaxy.quenching])
thubble_end = np.array([galaxy.t[1][i.below11] for i in galaxy.quenching])
quenching_times = np.array([i.quench_time for i in galaxy.quenching])
print(thubble_end,thubble_start)
print(galaxy.t[0])
print(galaxy.t[1])
print(galaxy.t[0][galaxy.quenching[0].indx])
print(galaxy.quenching[0].above9,galaxy.quenching[0].below11)
print(galaxy.quenching)
props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
for i in range(0, len(thubble_start)):
	axes[0].plot([float(thubble_start[i]),float(thubble_start[i])],[-12,-8], linestyle=':', color='b')
	axes[0].plot([float(thubble_end[i]),float(thubble_end[i])],[-12,-8], linestyle=':', color='r')
	xpos = thubble_start[i]-0.6
	axes[0].text(xpos, -9, r'$\tau_{q} = $'+'{:.3}'.format(quenching_times[i])+r' Gyr', fontsize=8, bbox=props)
for i in range(0, len(mergers_idx)):
	axes[0].plot(galaxy.t[0][mergers_idx[i]], np.log10(galaxy.ssfr[0][mergers_idx[i]]), marker='o', alpha=0.5, color='r', markersize=10)
for i in range(0, len(reju_id)):
	axes[0].plot(galaxy.t[0][reju_id[i]], np.log10(galaxy.ssfr[0][reju_id[i]]), marker='o', alpha=0.5, color='g', markersize=10)
axes[0].set_xlim([galaxy.t[0].min(),galaxy.t[0].max()])
axes[0].set_ylim([-11.5,-8])
axes[0].set_ylabel(r'$\log$(sSFR[yr$^{-1}$])', fontsize=16)
axes[0].legend(loc=1)


axes[1].plot(galaxy.t[0], np.log10(galaxy.m[0]), 'k-')
props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
for i in range(0, len(thubble_start)):
	axes[1].plot([thubble_start[i],thubble_start[i]],[np.log10(galaxy.m[0]).min(),np.log10(galaxy.m[0]).max()], linestyle=':', color='b')
	axes[1].plot([thubble_end[i],thubble_end[i]],[np.log10(galaxy.m[0]).min(),np.log10(galaxy.m[0]).max()], linestyle=':', color='r')
	xpos = thubble_start[i]-0.6
	axes[1].text(xpos, 10.5, r'$\tau_{q} = $'+'{:.3}'.format(quenching_times[i])+r' Gyr', fontsize=8, bbox=props)
for i in range(0, len(mergers_idx)):
	axes[1].plot(galaxy.t[0][mergers_idx[i]], np.log10(galaxy.m[0][mergers_idx[i]]), marker='o', alpha=0.5, color='r', markersize=10)
for i in range(0, len(reju_id)):
	axes[0].plot(galaxy.t[0][reju_id[i]], np.log10(galaxy.m[0][reju_id[i]]), marker='o', alpha=0.5, color='g', markersize=10)
axes[1].set_xlabel(r't (Gyr)', fontsize=16)
axes[1].set_ylabel(r'$\log(M_*[M_{\odot}]$)', fontsize=16)


cosmo = FlatLambdaCDM(H0=100*0.68, Om0=0.3, Ob0=0.04799952624117699,Tcmb0=2.73)  # set our cosmological parameters
axZ = axes[0].twiny()
axZ.set_xlim([galaxy.t[0].min(),galaxy.t[0].max()])
topticks1 = np.array([2,1,0.5,0.25,0])  # desired redshift labels
topticks2 = cosmo.age(topticks1).value  # tick locations in time
axZ.set_xticklabels(topticks1)
axZ.set_xticks(topticks2)
axZ.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
axZ.xaxis.set_label_position('top') # set the position of the second x-axis to top
axZ.set_xlabel('z', fontsize=16)
axZ.tick_params(labelsize=12)
for i in range(0, 2):
	axes[i].tick_params(labelsize=12)
fig.subplots_adjust(hspace=0)
fig.savefig('sfh_mh_'+str(GALAXY)+'.png', dpi=250, bbox_inches='tight')
