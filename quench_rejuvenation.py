#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 18:56:07 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import sys
simfolder = '../progen_analysis/m50n512'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
counterfile = 'galaxy_count.txt'#input('Text file with total number of galaxies per snapshot: ')
simname = 'm50n512'#input('SIMBA simulation version: ')
results_folder = '../quench_analysis/'+str(simname)+'/'
from import_progen import importApp
#from mergerFinder import myrunningmedian
from quenchingFinder import GalaxyData, quenchingFinder2, rejuvenation_rate_calculator, quenching_histogram

# Extract progen data from txt files
d, ngal = importApp(str(simfolder))
print('Total number of galaxies at z = 0: '+str(ngal))

#Store the galaxies sorted in objects of type GalaxyData
galaxies = []
for i in range(ngal):
    ssfr_gal = (d['sfr_gal' + str(i)]/d['m_gal'+str(i)])[::-1]+1e-14
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['galaxy_t' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['fgas_gal'+str(i)][::-1]
    gal_type = d['gal_type'+str(i)][::-1]
    galaxy = GalaxyData(i, ssfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type)
    galaxies.append(galaxy)

max_ngal = len(galaxies)
mass_limit = 9.5
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
        reju_m.append(np.log10(galaxy.rate[k+2]))

rates, red_cent, rates_sig = rejuvenation_rate_calculator(d, reju_z, counterfile)

print('Total number of rejuvenations: '+str(len(reju_z)))
print('Number of quenching events in first loop: '
        +str(sum([1 for galaxy in galaxies[0:max_ngal] for quench in galaxy.quenching])))
#Interpolation analysis
print("Interpolating")

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

scatter_labels = [['Final quenching Sat', 'Non-final quenching Sat', 'Final quenching Sat with rejuvenation' ],['Final quenching Central', 'Non-final quenching Central', 'Final quenching Central with rejuvenation']]
scatter_markers = ['.','*', '.']
scatter_faces = [['#e90c35','#ec5b75','none'],['#0c86e4','#8d0ce4','none']]
scatter_colours = [['#e90c35','#ec5b75','#e90c35'],['#0c86e4','#8d0ce4','#0c86e4']]
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
            ste_mass2_all.append(np.log10(galaxy.m_gal))
            quenching_times2_all.append(np.log10(quench.quench_time/galaxy.galaxy_t[end]))
            frac_gas2_all.append(np.log10(galaxy.fgas_gal))
            thubble2_all.append(np.log10(galaxy.galaxy_t[end]))
print(len(quenching_times2[0][0]), len(quenching_times2[0][1]), len(quenching_times2[0][2]), len(quenching_times2[0][0])+len(quenching_times2[0][1]))
print(len(quenching_times2[1][0]), len(quenching_times2[1][1]), len(quenching_times2[1][2]), len(quenching_times2[1][0])+len(quenching_times2[1][1]))
#print(galaxies_interpolated[0].quenching[-1])
print('Number of quenching events in second loop: '
        +str(sum([1 for galaxy in galaxies_interpolated for quench in galaxy.quenching])))

# tau = np.linspace(0.01,0.1, 100)
# sfr_tau = np.zeros(len(tau))
# for tau_i in range(0, len(tau)):
#     sfr_tau[tau_i] = 1/(80*0.0003*tau[tau_i])
# Plot the results
import seaborn as sns
sns.set(style="white")
fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
figre = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
# fig4 = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
# fig5 = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
# ax4 = fig4.add_subplot(1,1,1)
# ax5 = fig5.add_subplot(1,1,1)
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
axre = figre.add_subplot(1,1,1)
# ax5.plot(np.log10(tau), np.log10(sfr_tau), 'k-')
for i in range(0, len(scatter_labels)):
    for j in range(0, len(scatter_markers)):
        if j==2:
            ax1.scatter(redshifts2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            ax2.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            axre.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
        else:
            ax1.scatter(redshifts2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            ax2.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            axre.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
        # if scatter_labels[i][j]=='Non-final quenching Sat' or scatter_labels[i][j]=='Non-final quenching Central':
        #     print(quenching_times2[i][j])
        #print(len(quenching_times2[i][j]))
        # red = np.asarray(redshifts2[i][j])

        #ax4.scatter(np.log10(1+red), ste_mass2[i][j], c=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j])
        #ax5.scatter(np.log10(quenching_times2[i][j]), sfr_2[i][j], edgecolor=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j])
# ax1.plot((0.0, 4.0),(-1.5, -1.5), 'k--')
# ax2.plot((9.5, 11.5),(-1.5, -1.5), 'k--')
# ax2.plot((10.0, 10.0), (-3.2, 0.0), 'k--')

#x,y,ysig = myrunningmedian(red_cent, rates, 20)
#ax3.plot(x,y, 'k--')

ax1.set_xlabel(r'$z$', fontsize=16)
ax1.set_ylabel(r'$\log(t_{q})$', fontsize=16)
ax1.legend(loc='best', prop={'size': 10})
ax2.set_xlabel(r'$log(M_{*})$', fontsize=16)
ax2.set_ylabel(r'$\log(t_{q})$', fontsize=16)
fig.tight_layout()
fig.savefig(str(results_folder)+'quenching_morethan9.5_20.png', format='png', dpi=250)
axre.set_xlabel(r'$log(M_{*})$', fontsize=16)
axre.set_ylabel(r'$\log(t_{q}/t_{U})$', fontsize=16)
axre.legend(loc='best', prop={'size': 10})
figre.tight_layout()
figre.savefig(str(results_folder)+'quenching_morethan9.5_forreport.eps', format='eps', dpi=250)
# ax4.set_ylabel(r'$log(M_{*})$', fontsize=16)
# ax4.set_xlabel(r'$\log(1+z)$', fontsize=16)
# ax4.legend(loc='best', prop={'size': 10})
# fig4.tight_layout()
# fig4.savefig('quenching_morethan9.5_massred.png', format='png', dpi=250)
# ax5.set_ylabel(r'$log$(SFR (M$\odot$yr$^{-1}$))', fontsize=16)
# ax5.set_xlabel(r'$\log(t_{q})$', fontsize=16)
# ax5.legend(loc='best', prop={'size': 10})
# fig5.tight_layout()
# fig5.savefig('quenching_morethan9.5_sfrtime.png', format='png', dpi=250)

# fig = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
# ax3 = fig.add_subplot(1,1,1)
# ax3.set_xlabel(r'z', fontsize=16)
# ax3.set_ylabel(r'$\Gamma_{Rej}$ (Gyr$^{-1}$)', fontsize=16)
# leny = len(rates)-7
# ax3.errorbar(red_cent[0:leny], rates[0:leny], yerr=rates_sig[0:leny], capsize=2, marker='d', markerfacecolor='None', color='k', linestyle='--', markersize=12)
# fig.tight_layout()
# fig.savefig('rejuvenation_rate.png', format='png', dpi=250)

# fig2 = plt.figure(num=None, figsize=(10, 6), dpi=250, facecolor='w', edgecolor='k')
# fig3 = plt.figure(num=None, figsize=(10, 6), dpi=250, facecolor='w', edgecolor='k')
# quenchingPerType = {}
# for i in range(0, 3):
#     quenchingPerType['times'+str(i)] = []
#     quenchingPerType['redshifts'+str(i)] = []
#     quenchingPerType['mass'+str(i)] = []
#     quenchingPerType['f_gas'+str(i)] = []
#
# for j in range(0, len(quenching_times2_all)):
#     m = ste_mass2_all[j]
#     if 9.5<=m<10.3:
#         type = 0
#     elif 10.3<=m<11.0:
#         type = 1
#     elif m>=11.0:
#         type = 2
#     quenchingPerType['times'+str(type)].append(quenching_times2_all[j])
#     quenchingPerType['redshifts'+str(type)].append(redshifts2_all[j])
#     quenchingPerType['mass'+str(type)].append(m)
#     quenchingPerType['f_gas'+str(type)].append(frac_gas2_all[j])
# ax1 = fig2.add_subplot(1,1,1)
# ax2 = fig3.add_subplot(1,1,1)
# edgcolors = ['b', 'r', 'g']
# mass_ranges = [9.5,10.3,11.0,18.0]
# types = [r'$9.5\leq \log(M_*) < 10.3$', r'$10.3\leq \log(M_*) < 11.0$', r'$\log(M_*) \geq 11.0$']
# for i in range(0, len(mass_ranges)-1):
#     red_cent,frequency,frequency_sig,times,times_sig = quenching_histogram(d['z_gal0'],galaxies,max_ngal,mass_ranges[i],mass_ranges[i+1],quenchingPerType['times'+str(i)],
#                                                                             quenchingPerType['redshifts'+str(i)], 5)
#     ax1.errorbar(red_cent, frequency,yerr=frequency_sig, label=types[i], marker='o', linestyle='--', capsize=3, markersize=8)
#     ax2.errorbar(red_cent, times, yerr = times_sig, label=types[i], marker='o', linestyle='--', capsize=3, markersize=8)
# ax1.set_ylabel('Number of quenching events per galaxy', fontsize=16)
# ax1.set_xlabel('z', fontsize=16)
# ax2.set_ylabel(r' $\langle \log(t_{q}/t_{U}) \rangle$', fontsize=16)
# ax2.set_xlabel('z', fontsize=16)
# ax2.legend(loc='best', prop={'size': 12})
# ax1.legend(loc='best', prop={'size': 12})
# fig2.tight_layout()
# fig3.tight_layout()
# fig2.savefig('quenching_histogram_2.png', format='png', dpi=200)
# fig3.savefig('quenching_timeshisto.png', format='png', dpi=200)
#
# print('Quenching and Rejuvenation analysis done.')