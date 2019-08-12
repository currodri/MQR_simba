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
from astropy.cosmology import FlatLambdaCDM
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

# Import other codes
from mergerFinder import myrunningmedian, merger_finder
from quenchingFinder import quenchingFinder
from galaxy_class import GalaxyData, Quench
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
sfr_gal = d['sfr' + str(ourgalaxy_n)][::-1]
z_gal = d['z' + str(ourgalaxy_n)][::-1]
galaxy_t = d['t' + str(ourgalaxy_n)][::-1]
galaxy_m = d['m'+str(ourgalaxy_n)][::-1]
h2_gas = d['h2_gas'+str(ourgalaxy_n)][::-1]
h1_gas = 0
bh_m = 0
bhar = 0
local_den = 0
gal_type = d['g_type'+str(ourgalaxy_n)][::-1]
gal_pos = d['pos'+str(ourgalaxy_n)][::-1]
caesar_id = d['caesar_id'+str(ourgalaxy_n)][::-1]
galaxy = GalaxyData(ourgalaxy_n, sfr_gal, galaxy_m, z_gal, galaxy_t, h1_gas, h2_gas, bh_m, bhar,local_den, gal_type, gal_pos, caesar_id)
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
galaxies_interpolated = quenchingFinder(galaxies[0:max_ngal], 1, mass_limit)


print('Number of quenching events in first loop: '
		+str(sum([1 for galaxy in galaxies[0:max_ngal] for quench in galaxy.quenching])))
#Interpolation analysis
print("Interpolating...")

quenchingFinder(galaxies_interpolated, 1, mass_limit, True)

# Save results of rejuvenations coming from first loop
reju_z = []
reju_m = []
reju_t = []


for i in range(len(galaxies_interpolated)):
    galaxy = galaxies_interpolated[i]
    for index in galaxy.rejuvenations:
        reju_z.append(galaxy.z[index])
        reju_t.append(galaxy.t[0][index])
        reju_m.append(np.log10(galaxy.m[0][index]))

# Save quenchings in the classification schemed chosen
redshifts2_all = []
quenching_times2_all = []
ste_mass2_all = []
thubble2_all = []

redshifts2 = [[[],[],[]],[[],[],[]]]
quenching_times2 = [[[],[],[]],[[],[],[]]]
ste_mass2 = [[[],[],[]],[[],[],[]]]
thubble2 = [[[],[],[]],[[],[],[]]]


finalis = 0
nofinalis = 0
for i in range(0, len(galaxies_interpolated)):
    galaxy = galaxies_interpolated[i]
    if not isinstance(galaxy.m[1], int):
        if len(galaxy.quenching)>0:
            lastquench = galaxy.quenching[-1]
            # if len(galaxy.quenching)>1:
            #     for caca in range(0, len(galaxy.quenching)):
            #         print(galaxy.quenching[caca].quench_time)
        for quench in galaxy.quenching:
            start = quench.above9 + 1
            end = quench.below11
            if quench is lastquench:
                pos = 0
                finalis = finalis + 1
                for index in galaxy.rejuvenations:
                    if 0.2*galaxy.t[1][start]> (galaxy.t[1][start] - galaxy.t[0][index]) >=0:
                        pos = 2
            else:
                nofinalis = nofinalis + 1
                pos = 1
            #print(len(galaxy.z_gal), start, end, galaxy.id)
            if np.log10(galaxy.m[1][end])>=mass_limit:
                q_indx = int(quench.indx)
                q_type = 1 - int(galaxy.g_type[q_indx])
                redshifts2[q_type][pos].append(galaxy.z[q_indx])
                ste_mass2[q_type][pos].append(np.log10(galaxy.m[1][end]))
                quenching_times2[q_type][pos].append(np.log10(quench.quench_time/galaxy.t[1][end]))
                thubble2[q_type][pos].append(np.log10(galaxy.t[1][end]))
                redshifts2_all.append(galaxy.z[q_indx])
                ste_mass2_all.append(galaxy.m[1][end])
                quenching_times2_all.append(quench.quench_time)
                thubble2_all.append(galaxy.t[1][end])
print(len(quenching_times2[0][0]), len(quenching_times2[0][1]), len(quenching_times2[0][2]), len(quenching_times2[0][0])+len(quenching_times2[0][1]))
print(len(quenching_times2[1][0]), len(quenching_times2[1][1]), len(quenching_times2[1][2]), len(quenching_times2[1][0])+len(quenching_times2[1][1]))
print('Quenching and Rejuvenation analysis done.')
print(' ')
# Plot the results
import seaborn as sns
sns.set(style="white")
fig, axes = plt.subplots(2, 1, sharex='col', figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
axes[0].plot(galaxy_t, np.log10(ssfr_gal), 'k-')
axes[0].plot(galaxy_t, above, 'b--', label=r'Star-forming threshold: sSFR $=1/t_{U}$')
axes[0].plot(galaxy_t, below, 'r--', label=r'Quench threshold: sSFR $=0.2/t_{U}$')
axes[0].plot(galaxies_interpolated[0].galaxy_t, np.log10(galaxies_interpolated[0].ssfr_gal), linestyle='--', color='grey', alpha=0.7)
mergers_idx = np.asarray([np.where(galaxy_t==merg.galaxy_t[1])[0][0] for merg in mergers])
props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
axes[0].plot([8.739101250191442,8.739101250191442],[-12,-8], linestyle='-', color='k')
axes[0].plot([8.421918404720678,8.421918404720678],[-12,-8], linestyle='-.', color='k')
for i in range(0, len(thubble_start)):
	axes[0].plot([thubble_start[i],thubble_start[i]],[-12,-8], linestyle=':', color='b')
	axes[0].plot([thubble_end[i],thubble_end[i]],[-12,-8], linestyle=':', color='r')
	xpos = thubble_start[i]-0.6
	axes[0].text(xpos, -9, r'$\tau_{q} = $'+'{:.3}'.format(quenching_times[i])+r' Gyr', fontsize=8, bbox=props)
for i in range(0, len(mergers_idx)):
	axes[0].plot(mergers[i].galaxy_t[1], np.log10(ssfr_gal[mergers_idx[i]]), marker='o', alpha=0.5, color='r', markersize=10)
for i in range(0, len(reju_id)):
	axes[0].plot(reju_t[i], np.log10(galaxies_interpolated[0].ssfr_gal[reju_id[i]]), marker='o', alpha=0.5, color='g', markersize=10)
axes[0].set_xlim([galaxy_t.min(),galaxy_t.max()])
axes[0].set_ylim([-11.5,-8])
axes[0].set_ylabel(r'$\log$(sSFR[yr$^{-1}$])', fontsize=16)
axes[0].legend(loc=1)


axes[1].plot(galaxy_t, np.log10(galaxy_m), 'k-')
mergers_idx = np.asarray([np.where(galaxy_t==merg.galaxy_t[1])[0][0] for merg in mergers])
props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
for i in range(0, len(thubble_start)):
	axes[1].plot([thubble_start[i],thubble_start[i]],[np.log10(galaxy_m).min(),np.log10(galaxy_m).max()], linestyle=':', color='b')
	axes[1].plot([thubble_end[i],thubble_end[i]],[np.log10(galaxy_m).min(),np.log10(galaxy_m).max()], linestyle=':', color='r')
	xpos = thubble_start[i]-0.6
	axes[1].text(xpos, 10.5, r'$\tau_{q} = $'+'{:.3}'.format(quenching_times[i])+r' Gyr', fontsize=8, bbox=props)
for i in range(0, len(mergers_idx)):
	print(mergers[i].galaxy_t[2],np.log10(mergers[i].m_gal[2]))
	axes[1].plot(mergers[i].galaxy_t[2], np.log10(mergers[i].m_gal[2]), marker='o', alpha=0.5, color='r', markersize=10)
for i in range(0, len(reju_id)):
	axes[0].plot(reju_t[i], np.log10(galaxies_interpolated[0].m_gal[reju_id[i]]), marker='o', alpha=0.5, color='g', markersize=10)
axes[1].set_xlabel(r't (Gyr)', fontsize=16)
axes[1].set_ylabel(r'$\log(M_*[M_{\odot}]$)', fontsize=16)


cosmo = FlatLambdaCDM(H0=100*0.68, Om0=0.3, Ob0=0.04799952624117699,Tcmb0=2.73)  # set our cosmological parameters
axZ = axes[0].twiny()
axZ.set_xlim([galaxy_t.min(),galaxy_t.max()])
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
