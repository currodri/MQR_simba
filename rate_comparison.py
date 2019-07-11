#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 20 June 2019

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
from import_progen import importApp
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
for i in range(ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['t_gal' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]
    gal_type = d['gal_type'+str(i)][::-1]
    gal_pos = d['gal_pos'+str(i)][::-1]
    galaxy = GalaxyData(i, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type, gal_pos)
    galaxies.append(galaxy)

max_ngal = len(galaxies)
mass_limit = 9.5
min_merger_ratio = 0.2
max_redshift_mergers = 2.5

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

def Mass_Bin_Type(mass_bins, m_gal):
    mass_type = 3
    for mbin in range(0, len(mass_bins)):
        if 10**mass_bins[mbin][0] <= m_gal < 10**mass_bins[mbin][1]:
            mass_type = mbin
    return mass_type

def Fractional_Rate(mergers,sf_galaxies,q_masses,q_reds,q_thubble,reju_z,reju_t,reju_m,n_bins,max_redshift_mergers):
    mass_limits = [[9.5,10.3], [10.3,11.0],[11.0,18.0]]
    mass_labels = [r'$9.5\leq \log(M_*) < 10.3$', r'$10.3\leq \log(M_*) < 11.0$', r'$\log(M_*) \geq 11.0$']
    z_bins = np.linspace(0.0, max_redshift_mergers, n_bins)
    r_merger = {}
    r_quench = {}
    r_reju = {}
    r_merger['all'] = np.zeros(n_bins-1)
    r_quench['all'] = np.zeros(n_bins-1)
    r_reju['all'] = np.zeros(n_bins-1)
    for bini in range(0, len(mass_limits)):
        r_merger['massbin'+str(bini)] = np.zeros(n_bins-1)
        r_quench['massbin'+str(bini)] = np.zeros(n_bins-1)
        r_reju['massbin'+str(bini)] = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        sf_counter = {}
        for sf_type in range(0, len(mass_limits)):
            sf_counter['massbin'+str(sf_type)] = 0
        counter = 0
        times = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[2] < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits,merger.m_gal[2])
                if type != 3:
                    if type==0:
                        counter = counter + 1
                    r_merger['massbin'+str(type)][i] = r_merger['massbin'+str(type)][i] + 1
                    times.append(merger.galaxy_t[2])
        for k in range(0, len(sf_galaxies)):
            sf = sf_galaxies[k]
            if z_bins[i]<= sf.z_gal < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits,sf.m_gal)
                if type != 3:
                    sf_counter['massbin'+str(type)] = sf_counter['massbin'+str(type)] + 1
                    times.append(sf.galaxy_t)
        for l in range(0, len(q_reds)):
            quench_red = q_reds[l]
            if z_bins[i]<= quench_red < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits,q_masses[l])
                if type != 3:
                    r_quench['massbin'+str(type)][i] = r_quench['massbin'+str(type)][i] + 1
                    times.append(q_thubble[l])
        for m in range(0, len(reju_z)):
            reju = reju_z[m]
            if z_bins[i]<= reju < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits, reju_m[m])
                if type != 3:
                    r_reju['massbin'+str(type)][i] = r_reju['massbin'+str(type)][i] + 1
                    times.append(reju_t[m])
        times = np.asarray(times)
        delta_t = float(times.max() - times.min())
        a = 0
        b = 0
        c = 0
        sf = 0
        for ty in range(0, len(mass_limits)):
            #print(counter)
            #print(r_merger['massbin'+str(ty)][i], r_quench['massbin'+str(ty)][i], r_merger['massbin'+str(ty)][i])
            normalization = float(float(r_merger['massbin'+str(ty)][i]+sf_counter['massbin'+str(ty)])*delta_t)
            a = a + r_merger['massbin'+str(ty)][i]
            b = b + r_quench['massbin'+str(ty)][i]
            c = c + r_reju['massbin'+str(ty)][i]
            sf = sf + sf_counter['massbin'+str(ty)]
            r_merger['massbin'+str(ty)][i] = float(r_merger['massbin'+str(ty)][i])/normalization
            r_quench['massbin'+str(ty)][i] = float(r_quench['massbin'+str(ty)][i])/normalization
            r_reju['massbin'+str(ty)][i] = float(r_reju['massbin'+str(ty)][i])/normalization
        normalization = float(float(a+sf)*delta_t)
        r_merger['all'][i] = float(a)/normalization
        r_quench['all'][i] = float(b)/normalization
        r_reju['all'][i] = float(c)/normalization
    fig, ax = plt.subplots(3, 1, sharex='col', num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    x_dat = np.log10(1+z_cent)
    for i in range(0, len(mass_limits)):
        ax[0].plot(x_dat, np.log10(r_merger['massbin'+str(i)]), linestyle='--', marker='d', label=mass_labels[i])
        ax[1].plot(x_dat, np.log10(r_quench['massbin'+str(i)]), linestyle='--', marker='d')
        ax[2].plot(x_dat, np.log10(r_reju['massbin'+str(i)]), linestyle='--', marker='d')

    ax[0].plot(x_dat, np.log10(r_merger['all']),linestyle='--', marker='d', color='k', label='All')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_dat, np.log10(r_merger['all']))
    ax[0].plot(x_dat,np.log10((10**intercept)*(1+z_cent)**(slope)), 'k-', label=r'$10^{%.2f}\cdot(1+z)^{%.2f}$' % (intercept, slope) )
    print("slope: %f    intercept: %f    r_value: %f    p_value: %f    std_error: %f" % (slope, intercept,r_value, p_value, std_err))

    ax[1].plot(x_dat, np.log10(r_quench['all']),linestyle='--', marker='d', color='k')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_dat, np.log10(r_quench['all']))
    ax[1].plot(x_dat,np.log10((10**intercept)*(1+z_cent)**(slope)), 'k-', label=r'$10^{%.2f}\cdot(1+z)^{%.2f}$' % (intercept, slope) )
    print("slope: %.f    intercept: %f    r_value: %f    p_value: %f    std_error: %f" % (slope, intercept,r_value, p_value, std_err))

    ax[2].plot(x_dat, np.log10(r_reju['all']),linestyle='--', marker='d', color='k')
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_dat, np.log10(r_reju['all']))
    ax[2].plot(x_dat,np.log10((10**intercept)*(1+z_cent)**(slope)), 'k-', label=r'$10^{%.2f}\cdot(1+z)^{%.2f}$' % (intercept, slope) )
    print("slope: %f    intercept: %f    r_value: %f    p_value: %f    std_error: %f" % (slope, intercept,r_value, p_value, std_err))
    ax[0].set_ylabel(r'$\log(\mathcal{R}_{Mer})$ [Gyr$^{-1}$]', fontsize=16)
    ax[1].set_ylabel(r'$\log(\mathcal{R}_{Que})$ [Gyr$^{-1}$]', fontsize=16)
    ax[2].set_ylabel(r'$\log(\mathcal{R}_{Rej})$ [Gyr$^{-1}$]', fontsize=16)
    ax[2].set_xlabel(r'$\log(1+z)$', fontsize=16)
    axR = ax[0].twiny()
    maxlz = 0.56
    ax[0].set_xlim(0.03,maxlz)
    axR.set_xlim(0.03,maxlz)
    topticks1 = np.array([0,1,2])  # desired redshift labels
    topticks2 = np.log10(1+topticks1)  # tick locations in time
    axR.set_xticklabels(topticks1)
    axR.set_xticks(topticks2)
    axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
    axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
    axR.set_xlabel('z', fontsize=16)
    axR.tick_params(labelsize=12)
    ax[0].legend(loc='best', prop={'size': 12}, fontsize=14)
    ax[1].legend(loc='best', prop={'size': 12}, fontsize=14)
    ax[2].legend(loc='best', prop={'size': 12}, fontsize=14)
    for i in range(0,3):
        ax[i].tick_params(labelsize=12)
    fig.subplots_adjust(hspace=0)
    fig.savefig(str(results_folder)+'mqr_fractional_rate.png', format='png', dpi=200, bbox_inches='tight')

def Density_Rate(mergers,q_masses,q_reds,q_thubble,reju_z,reju_t,reju_m,n_bins,max_redshift_mergers):
    mass_limits = [[9.5,10.3], [10.3,11.0],[11.0,18.0]]
    mass_labels = [r'$9.5\leq \log(M_*) < 10.3$', r'$10.3\leq \log(M_*) < 11.0$', r'$\log(M_*) \geq 11.0$']
    z_bins = np.linspace(0.0, max_redshift_mergers, n_bins)
    r_merger = {}
    r_quench = {}
    r_reju = {}
    for bini in range(0, len(mass_limits)):
        r_merger['massbin'+str(bini)] = np.zeros(n_bins-1)
        r_quench['massbin'+str(bini)] = np.zeros(n_bins-1)
        r_reju['massbin'+str(bini)] = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        counter = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits,merger.m_gal[1])
                if type != 3:
                    r_merger['massbin'+str(type)][i] = r_merger['massbin'+str(type)][i] + 1
                    counter.append(merger.z_gal[1])
        for l in range(0, len(q_reds)):
            quench_red = q_reds[l]
            if z_bins[i]<= quench_red < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits,q_masses[l])
                if type != 3:
                    r_quench['massbin'+str(type)][i] = r_quench['massbin'+str(type)][i] + 1
                    counter.append(quench_red)
        for m in range(0, len(reju_z)):
            reju = reju_z[m]
            if z_bins[i]<= reju < z_bins[i+1]:
                type = Mass_Bin_Type(mass_limits, reju_m[m])
                if type != 3:
                    r_reju['massbin'+str(type)][i] = r_reju['massbin'+str(type)][i] + 1
                    counter.append(reju)
        total_v = float(100*len(np.unique(counter)))
        for ty in range(0, len(mass_limits)):
            r_merger['massbin'+str(ty)][i] = float(r_merger['massbin'+str(ty)][i])/total_v
            r_quench['massbin'+str(ty)][i] = float(r_quench['massbin'+str(ty)][i])/total_v
            r_reju['massbin'+str(ty)][i] = float(r_reju['massbin'+str(ty)][i])/total_v
    fig, ax = plt.subplots(3, 1, sharex='col', num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    x_dat = np.log10(1+z_cent)
    for i in range(0, len(mass_limits)):
        ax[0].plot(x_dat, np.log10(r_merger['massbin'+str(i)]), linestyle='--', marker='d', label=mass_labels[i])
        ax[1].plot(x_dat, np.log10(r_quench['massbin'+str(i)]), linestyle='--', marker='d')
        ax[2].plot(x_dat, np.log10(r_reju['massbin'+str(i)]), linestyle='--', marker='d')
    ax[0].set_ylabel(r'$\log(\Gamma_{Mer})$', fontsize=16)
    ax[1].set_ylabel(r'$\log(\Gamma_{Que})$', fontsize=16)
    ax[2].set_ylabel(r'$\log(\Gamma_{Rej})$', fontsize=16)
    ax[2].set_xlabel(r'$\log(1+z)$', fontsize=16)
    axR = ax[0].twiny()
    maxlz = 0.56
    ax[0].set_xlim(0.03,maxlz)
    axR.set_xlim(0.03,maxlz)
    topticks1 = np.array([0,1,2])  # desired redshift labels
    topticks2 = np.log10(1+topticks1)  # tick locations in time
    axR.set_xticklabels(topticks1)
    axR.set_xticks(topticks2)
    axR.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
    axR.xaxis.set_label_position('top') # set the position of the second x-axis to top
    axR.set_xlabel('z', fontsize=16)
    axR.tick_params(labelsize=12)
    for i in range(0,3):
        ax[i].tick_params(labelsize=12)
    ax[0].legend(loc='best', prop={'size': 12}, fontsize=14)
    fig.subplots_adjust(hspace=0)
    fig.savefig(str(results_folder)+'mqr_density_rate.png', format='png', dpi=200, bbox_inches='tight')

Fractional_Rate(mergers,sf_galaxies,ste_mass2_all,redshifts2_all,thubble2_all,reju_z,reju_t,reju_m,10,max_redshift_mergers)
Density_Rate(mergers,ste_mass2_all,redshifts2_all,thubble2_all,reju_z,reju_t,reju_m,10,max_redshift_mergers)
