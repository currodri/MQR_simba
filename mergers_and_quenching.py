#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 June 2019

@author: currorodriguez
"""
# Import required libraries
import numpy as np
import matplotlib
from scipy import stats
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")

# Import other codes
from import_progen import importApp
from mergerFinder import myrunningmedian, merger_finder
from quenchingFinder import GalaxyData, quenchingFinder2, rejuvenation_rate_calculator, quenching_histogram
import sys
simfolder = '../progen_analysis/m100n1024'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
counterfile = '../progen_analysis/m100n1024/galaxy_count_m100n1024.txt'#input('Text file with total number of galaxies per snapshot: ')
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../mandq_relations/'+str(simname)+'/'
timefile = '../quench_analysis'+str(simname)+'/times_m100n1024.txt'
redshiftfile = '../quench_analysis'+str(simname)+'/redshifts_m100n1024.txt'

# Extract progen data from txt files
d, ngal = importApp(str(simfolder))
print('Total number of galaxies at z = 0: '+str(ngal))

#Store the galaxies sorted in objects of type GalaxyData
galaxies = []
for i in range(ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['galaxy_t' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]
    gal_type = d['gal_type'+str(i)][::-1]
    galaxy = GalaxyData(i, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type)
    galaxies.append(galaxy)

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

redshifts2_all = []
quenching_times2_all = []
ste_mass2_all = []
frac_gas2_all = []
thubble2_all = []

redshifts2 = [[[],[],[]],[[],[],[]]]
quenching_times2 = [[[],[],[]],[[],[],[]]]
ste_mass2 = [[[],[],[]],[[],[],[]]]
frac_gas2 = [[[],[],[]],[[],[],[]]]
thubble2 = [[[],[],[]],[[],[],[]]]
sfr_2 = [[[],[],[]],[[],[],[]]]


finalis = 0
nofinalis = 0
for i in range(0, len(galaxies_interpolated)):
    galaxy = galaxies_interpolated[i]
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
            for k in range(0, len(galaxy.rate), 3):
                if 0.2*galaxy.galaxy_t[start]> (galaxy.galaxy_t[start] - galaxy.rate[k+1]) >=0:
                    pos = 2
        else:
            nofinalis = nofinalis + 1
            pos = 1
        #print(len(galaxy.z_gal), start, end, galaxy.id)
        if np.log10(galaxy.m_gal)>=mass_limit:
            redshifts2[int(quench.type)][pos].append(galaxy.z_gal)
            ste_mass2[int(quench.type)][pos].append(np.log10(galaxy.m_gal))
            quenching_times2[int(quench.type)][pos].append(np.log10(quench.quench_time))
            frac_gas2[int(quench.type)][pos].append(np.log10(galaxy.fgas_gal))
            thubble2[int(quench.type)][pos].append(np.log10(galaxy.galaxy_t[end]))
            sfr_2[int(quench.type)][pos].append(np.log10(galaxy.ssfr_gal[start-1]*galaxy.m_gal))
            redshifts2_all.append(galaxy.z_gal)
            ste_mass2_all.append(galaxy.m_gal)
            quenching_times2_all.append(quench.quench_time)
            frac_gas2_all.append(galaxy.fgas_gal)
            thubble2_all.append(galaxy.galaxy_t[end])
print(len(quenching_times2[0][0]), len(quenching_times2[0][1]), len(quenching_times2[0][2]), len(quenching_times2[0][0])+len(quenching_times2[0][1]))
print(len(quenching_times2[1][0]), len(quenching_times2[1][1]), len(quenching_times2[1][2]), len(quenching_times2[1][0])+len(quenching_times2[1][1]))
print('Number of quenching events in second loop: '
        +str(sum([1 for galaxy in galaxies_interpolated for quench in galaxy.quenching])))
print('Quenching and Rejuvenation analysis done.')
print(' ')


def mergerquench_relation():
    print('Start finding for connection between mergers and quenching')
    time_diff = []
    merger_ratios = []
    quenching_times = []
    for i in range(0, len(mergers)):
        merg = mergers[i]
        possible_q = []
        possible_q_time = []
        for j in range(0, len(galaxies_interpolated)):
            galaxy = galaxies_interpolated[j]
            if galaxy.id == merg.id:
                for quench in galaxy.quenching:
                    start = quench.above9 #+ 1
                    end = quench.below11
                    if np.log10(galaxy.m_gal)>=mass_limit:
                        possible_q.append(galaxy.galaxy_t[start])
                        possible_q_time.append(quench.quench_time/galaxy.galaxy_t[end])
        diff = []
        for k in range(0, len(possible_q)):
            diff.append(possible_q[k] - merg.galaxy_t[1])
        diff = np.asarray(diff)
        #print(diff)
        if len(possible_q)>0:
            if diff[np.argmin(diff)]>=0:
                merger_ratios.append(merg.merger_ratio)
                time_diff.append(possible_q[np.argmin(diff)]-merg.galaxy_t[1])
                quenching_times.append(possible_q_time[np.argmin(diff)])
    time_diff = np.asarray(time_diff)
    merger_ratios = np.asarray(merger_ratios)
    quenching_times = np.asarray(quenching_times)
    fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r'$T_q - T_m$(Gyr)', fontsize=16)
    ax.set_ylabel(r'$\log(N) $(Gyr)', fontsize=16)
    ax.hist(time_diff, bins=12, histtype='step', log=True)
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mergertime_and_quenching.png',format='png', dpi=250)
    return time_diff,quenching_times,merger_ratios

def quench_delay(delay,quenching_times, merger_ratios):
    print('Making scatter plot of quenching times vs delay of the quenching process after a merger...')
    fig1 = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig1.add_subplot(1,1,1)
    pear = stats.pearsonr(delay, np.log10(quenching_times))
    print('The Pearson correlation R is given by: '+str(pear))
    sc1 = ax1.scatter(delay, np.log10(quenching_times), c=merger_ratios , cmap='viridis')
    ax1.set_xlabel(r'$T_q - T_m$', fontsize=16)
    ax1.set_ylabel(r'$\log(t_{quench}/t_{Hubble})$', fontsize=16)
    #ax1.set_xlim([0.18, 0.6])
    fig1.colorbar(sc1, ax=ax1, label=r'Merger ratio $R$', orientation='horizontal')
    fig1.tight_layout()
    fig1.savefig(str(results_folder)+'delay_and_qtime.png',format='png', dpi=250)

def merger_reju_relation():
    print('Start finding for connection between mergers and rejuvenations')
    time_diff = []
    merger_boost = []

    for i in range(0, len(mergers)):
        merg = mergers[i]
        possible_r = []
        for j in range(0, len(reju_t)):
            if reju_id[j]==merg.id:
                possible_r.append(reju_t[j])
        diff = []
        for k in range(0, len(possible_r)):
            diff.append(possible_r[k] - merg.galaxy_t[1])
        diff = np.asarray(diff)
        if len(possible_r)>0:
            if diff[np.argmin(diff)]>=0:
                if merg.fgas_boost<0:
                    merger_boost.append(0.001)
                else:
                    merger_boost.append(merg.fgas_boost)
                if merg.fgas_boost>10.0:
                    print(merg.id)
                time_diff.append(possible_r[np.argmin(diff)]-merg.galaxy_t[1])
    time_diff = np.asarray(time_diff)
    merger_boost = np.asarray(merger_boost)
    fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r'$T_r - T_m$(Gyr)', fontsize=16)
    ax.set_ylabel(r'$\log(N) $(Gyr)', fontsize=16)
    ax.hist(time_diff, bins=12, histtype='step', log=True)
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mergertime_and_rejuvenation.png',format='png', dpi=250)

time_diff, q_times, m_ratios = mergerquench_relation()
quench_delay(time_diff,q_times,m_ratios)
merger_reju_relation()
