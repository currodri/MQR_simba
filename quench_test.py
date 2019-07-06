#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 6 18:56:07 2019

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import matplotlib
import pickle
from scipy import stats
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

# Import other codes
from mergerFinder import myrunningmedian, merger_finder
from quenchingFinder import GalaxyData, quenchingFinder2, rejuvenation_rate_calculator, quenching_histogram
import sys
simfolder = '../progen_analysis/m100n1024'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
counterfile = '../progen_analysis/m100n1024/galaxy_count_m100n1024.txt'#input('Text file with total number of galaxies per snapshot: ')
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../rate_analysis/'+str(simname)+'/'
timefile = '../quench_analysis'+str(simname)+'/times_m100n1024.txt'
redshiftfile = '../quench_analysis'+str(simname)+'/redshifts_m100n1024.txt'
pickle_file = '../progen_analysis/m100n1024/progen_'+str(simname)+'.pkl'

# Extract progen data from txt files
# d, ngal = importApp(str(simfolder))
# print('Total number of galaxies at z = 0: '+str(ngal))
obj = open(pickle_file, 'rb')
d = pickle.load(obj)
ngal = 49215
print('Total number of galaxies at z = 0: '+str(ngal))

#Store the galaxies sorted in objects of type GalaxyData
galaxies = []
ourgalaxy_n = int(input('What galaxy do you wanna plot: '))
sfr_gal = d['sfr_gal' + str(ourgalaxy_n)][::-1]
ssfr_gal = (d['sfr_gal' + str(ourgalaxy_n)]/d['m_gal'+str(ourgalaxy_n)])[::-1]+1e-14
sfe_gal = d['sfe_gal' + str(ourgalaxy_n)][::-1]
z_gal = d['z_gal' + str(ourgalaxy_n)][::-1]
galaxy_t = d['t_gal' + str(ourgalaxy_n)][::-1]
galaxy_m = d['m_gal'+str(ourgalaxy_n)][::-1]
fgas_gal = d['h2_gal'+str(ourgalaxy_n)][::-1]
gal_type = d['gal_type'+str(ourgalaxy_n)][::-1]
gal_pos = d['gal_pos'+str(ourgalaxy_n)][::-1]
galaxy = GalaxyData(ourgalaxy_n, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type, gal_pos)
galaxies.append(galaxy)

above = []
below = []
for i in range(0, len(galaxy_t)):
    above.append(np.log10(1/galaxy_t[i])-9)
    below.append(np.log10(0.2/galaxy_t[i])-9)

max_ngal = len(galaxies)
mass_limit = 9.5
min_merger_ratio = 0.2
max_redshift_mergers = 3.5

# Perform the search for mergers
mergers, sf_galaxies = merger_finder(galaxies, min_merger_ratio, 10**mass_limit, max_redshift_mergers)

print('Merger analysis done.')

# Perform the quenching and rejuvenation analysis
galaxies_interpolated = quenchingFinder2(galaxies[0:max_ngal], 1, mass_limit)

# Save results of rejuvenations coming from first loop
reju_z = []
reju_m = []
reju_t = []


for i in range(len(galaxies)):
    galaxy = galaxies[i]
    for k in range(0, len(galaxy.rate), 3):
        #if np.log10(galaxy.rate[k+1])>=mass_limit:
        reju_z.append(galaxy.rate[k])
        reju_t.append(galaxy.rate[k+1])
        reju_m.append(galaxy.rate[k+2])


print('Total number of rejuvenations: '+str(len(reju_z)))
print('Number of quenching events in first loop: '
        +str(sum([1 for galaxy in galaxies[0:max_ngal] for quench in galaxy.quenching])))
#Interpolation analysis
print("Interpolating...")

quenchingFinder2(galaxies_interpolated, 1, mass_limit, True)

redshifts = []
ste_mass = []
quenching_times = []
frac_gas = []
thubble_start = []
thubble_end = []

for i in range(0, len(galaxies_interpolated)):
    galaxy = galaxies_interpolated[i]
    for quench in galaxy.quenching:
        start = quench.above9 + 1
        end = quench.below11
        redshifts.append(galaxy.z_gal)
        ste_mass.append(galaxy.m_gal)
        quenching_times.append(quench.quench_time)
        frac_gas.append(galaxy.fgas_gal)
        thubble_start.append(galaxy.galaxy_t[start])
        thubble_end.append(galaxy.galaxy_t[end])

print(redshifts, quenching_times, thubble_start, thubble_end)
print([merg.galaxy_t[1] for merg in mergers])
print(reju_t)
print('Number of quenching events in second loop: '
        +str(sum([1 for galaxy in galaxies_interpolated for quench in galaxy.quenching])))
# Plot the results
import seaborn as sns
sns.set(style="white")
fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
ax1 = fig.add_subplot(1,1,1)
ax1.plot(galaxy_t, np.log10(ssfr_gal), 'k-')
ax1.plot(galaxy_t, above, 'b--', label=r'Star-forming threshold: sSFR $=1/t_{U}$')
ax1.plot(galaxy_t, below, 'r--', label=r'Quench threshold: sSFR $=0.2/t_{U}$')
ax1.plot(galaxies_interpolated[0].galaxy_t, np.log10(galaxies_interpolated[0].ssfr_gal), 'r--')
mergers_idx = np.asarray([np.where(galaxy_t==merg.galaxy_t[1])[0][0] for merg in mergers])
rejuvenations_idx = np.asarray([np.where(galaxy_t==rej)[0][0] for rej in reju_t])
print(mergers_idx, rejuvenations_idx)
for i in range(0, len(thubble_start)):
    ax1.plot([thubble_start[i],thubble_start[i]],[-12,-8], linestyle=':', color='b')
    ax1.plot([thubble_end[i],thubble_end[i]],[-12,-8], linestyle=':', color='r')
for i in range(0, len(mergers_idx)):
    ax1.plot(mergers[i].galaxy_t[1], np.log10(ssfr_gal[mergers_idx[i]]), marker='o', alpha=0.5, color='r', markersize=10)
for i in range(0, len(rejuvenations_idx)):
    ax1.plot(reju_t[i], np.log10(ssfr_gal[rejuvenations_idx[i]]), marker='o', alpha=0.5, color='g', markersize=10)
# ax1.plot(7.668095312278215, np.log10(1.9816408261273213e-10), marker='o', alpha=0.5, color='g', markersize=10)
# ax1.text(2.0, -12.8, r'$t_{q} = $'+'{:.3}'.format(quenching_times[0])+r' Gyr')
# ax1.text(6.0, -12.8, r'$t_{q} = $'+'{:.3}'.format(quenching_times[1])+r' Gyr')
ax1.set_xlim([galaxy_t.min(),12])
ax1.set_ylim([-12,-8])
ax1.set_xlabel(r't (Gyr)', fontsize=16)
ax1.set_ylabel(r'$\log$ sSFR ($M_'+u'\u2609'+r'$yr$^{-1}$)', fontsize=16)
ax1.legend(loc=1)
fig.tight_layout()
fig.savefig('quench_finder_test.png', dpi=250)
