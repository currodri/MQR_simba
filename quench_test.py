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

# Extract quench data from pickle file
data_file = '/home/curro/quenchingSIMBA/code/SH_Project/mandq_results_%s.pkl' % (MODEL)
obj = open(data_file, 'rb')
d = pickle.load(obj)

galaxy = d['galaxies'][GALAXY]

above = np.zeros(len(galaxy.t[0]))
below = np.zeros(len(galaxy.t[0]))

for i in range(0, len(galaxy.t[0])):
	above[i] = np.log10(1/(galaxy.t[0][i]))-9
	below[i] = np.log10(0.2/(galaxy.t[0][i])) - 9


# Plot the results

fig, axes = plt.subplots(2, 1, sharex='col', figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
axes[0].plot(galaxy.t[0], np.log10(galaxy.ssfr[0]), 'k-')
axes[0].plot(galaxy.t[0], above, 'b--', label=r'Star-forming threshold: sSFR $=1/t_{U}$')
axes[0].plot(galaxy.t[0], below, 'r--', label=r'Quench threshold: sSFR $=0.2/t_{U}$')
axes[0].plot(galaxy.t[1], np.log10(galaxy.ssfr[1]), linestyle='--', color='grey', alpha=0.7)
mergers_idx = np.array([i.indx for i in galaxy.mergers])
reju_id = np.array([i for i in galaxy.rejuvenations])
thubble_start = np.array([galaxy.t[1][i.above9] for i in galaxy.quenching])
thubble_end = np.array([galaxy.t[1][i.below11] for i in galaxy.quenching])
quenching_times = np.array([i.quench_time for i in galaxy.quenching])
props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
axes[0].plot([8.739101250191442,8.739101250191442],[-12,-8], linestyle='-', color='k')
axes[0].plot([8.421918404720678,8.421918404720678],[-12,-8], linestyle='-.', color='k')
for i in range(0, len(thubble_start)):
	axes[0].plot([thubble_start[i],thubble_start[i]],[-12,-8], linestyle=':', color='b')
	axes[0].plot([thubble_end[i],thubble_end[i]],[-12,-8], linestyle=':', color='r')
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
	print(galaxy.t[0][mergers_idx[i]],np.log10(galaxy.m[0]))
	axes[1].plot(galaxy.t[0], np.log10(galaxy.m[0]), marker='o', alpha=0.5, color='r', markersize=10)
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
fig.savefig('quench_finder_test.png', dpi=250, bbox_inches='tight')
