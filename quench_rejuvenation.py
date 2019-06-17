#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 18:56:07 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
# Import other codes
from import_progen import importApp
from mergerFinder import myrunningmedian
from quenchingFinder import GalaxyData, quenchingFinder2, rejuvenation_rate_calculator, quenching_histogram
import sys
simfolder = '../progen_analysis/m100n1024'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
counterfile = '../progen_analysis/m100n1024/galaxy_count_m100n1024.txt'#input('Text file with total number of galaxies per snapshot: ')
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../quench_analysis/'+str(simname)+'/'
timefile = results_folder+'/times_m100n1024.txt'
redshiftfile = results_folder+'/redshifts_m100n1024.txt'


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
            ste_mass2_all.append(np.log10(galaxy.m_gal))
            quenching_times2_all.append(np.log10(quench.quench_time/galaxy.galaxy_t[end]))
            frac_gas2_all.append(np.log10(galaxy.fgas_gal))
            thubble2_all.append(np.log10(galaxy.galaxy_t[end]))
print(len(quenching_times2[0][0]), len(quenching_times2[0][1]), len(quenching_times2[0][2]), len(quenching_times2[0][0])+len(quenching_times2[0][1]))
print(len(quenching_times2[1][0]), len(quenching_times2[1][1]), len(quenching_times2[1][2]), len(quenching_times2[1][0])+len(quenching_times2[1][1]))
print('Number of quenching events in second loop: '
        +str(sum([1 for galaxy in galaxies_interpolated for quench in galaxy.quenching])))
print('Quenching and Rejuvenation analysis done.')
print(' ')
print('Now, you can choose what plots to obtain!')
print(' ')
print(' ')
# Plot the results

def Quenching_Scatter_Plot(redshifts2, quenching_times2, ste_mass2):
    scatter_labels = [['Final quenching Sat', 'Non-final quenching Sat', 'Final quenching Sat with rejuvenation' ],['Final quenching Central', 'Non-final quenching Central', 'Final quenching Central with rejuvenation']]
    scatter_markers = ['.','*', '.']
    cmaps = ['Reds', 'Blues']
    grids = [30,70]
    #scatter_colours = [['#e90c35','#ec5b75','#e90c35'],['#0c86e4','#8d0ce4','#0c86e4']]
    scatter_colours = [['#e90c35','m','r'],['#0c86e4','c','b']]
    fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    #figre = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
    # fig4 = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
    # fig5 = plt.figure(num=None, figsize=(7, 5), dpi=80, facecolor='w', edgecolor='k')
    # ax4 = fig4.add_subplot(1,1,1)
    # ax5 = fig5.add_subplot(1,1,1)
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #axre = figre.add_subplot(1,1,1)
    # ax5.plot(np.log10(tau), np.log10(sfr_tau), 'k-')
    for i in range(len(scatter_labels), 0, -1):
        for j in range(0, len(scatter_markers)):
            if j==0:
                hb2 = ax2.hexbin(ste_mass2[i][j], quenching_times2[i][j], cmap=cmaps[i], alpha=0.85, bins='log', gridsize=grids[i])
                hb1 = ax1.hexbin(redshifts2[i][j], quenching_times2[i][j], cmap=cmaps[i], alpha=0.75, bins='log', gridsize=grids[i])
                cb1 = fig.colorbar(hb1, ax=ax1)
                cb1.set_label(r'$\log$(N)')
                cb2 = fig.colorbar(hb2, ax=ax2)
                cb2.set_label(r'$\log$(N)')
            # if j==2:
            #     ax1.scatter(redshifts2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            #     ax2.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            #     #axre.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
            else:
                ax1.scatter(redshifts2[i][j], quenching_times2[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], edgecolor=scatter_colours[i][j], s=40, facecolors='w')
                ax2.scatter(ste_mass2[i][j], quenching_times2[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], edgecolor=scatter_colours[i][j], s=40, facecolors='w')
                #axre.scatter(ste_mass2[i][j], quenching_times2[i][j], edgecolor=scatter_colours[i][j], label =scatter_labels[i][j], marker=scatter_markers[j], facecolors=scatter_faces[i][j], s=60)
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
    #axre.set_xlabel(r'$log(M_{*})$', fontsize=16)
    #axre.set_ylabel(r'$\log(t_{q}/t_{U})$', fontsize=16)
    #axre.legend(loc='best', prop={'size': 10})
    #figre.tight_layout()
    #figre.savefig(str(results_folder)+'quenching_morethan9.5_forreport.eps', format='eps', dpi=250)
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

def Rejuvenation_Rate_Plot(d, reju_z, counterfile, timefile, redshiftfile):
    rates, red_cent, rates_sig = rejuvenation_rate_calculator(d, reju_z, counterfile, timefile, redshiftfile)
    fig = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax3 = fig.add_subplot(1,1,1)
    ax3.set_xlabel(r'z', fontsize=16)
    ax3.set_ylabel(r'$\Gamma_{Rej}$ (Gyr$^{-1}$)', fontsize=16)
    leny = len(rates)#-7
    ax3.errorbar(red_cent[0:leny], rates[0:leny], yerr=rates_sig[0:leny], capsize=2, marker='d', markerfacecolor='None', color='k', linestyle='--', markersize=12)
    fig.tight_layout()
    fig.savefig(str(results_folder)+'rejuvenation_rate.png', format='png', dpi=250)

def Quenching_Histogram_Plots(quenching_times2_all, redshifts2_all, ste_mass2_all, frac_gas2_all):
    fig2, axes = plt.subplots(2, 1, sharex='col', num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
    quenchingPerType = {}
    for i in range(0, 3):
        quenchingPerType['times'+str(i)] = []
        quenchingPerType['redshifts'+str(i)] = []
        quenchingPerType['mass'+str(i)] = []
        quenchingPerType['f_gas'+str(i)] = []

    for j in range(0, len(quenching_times2_all)):
        m = ste_mass2_all[j]
        if 9.5<=m<10.3:
            type = 0
        elif 10.3<=m<11.0:
            type = 1
        elif m>=11.0:
            type = 2
        quenchingPerType['times'+str(type)].append(quenching_times2_all[j])
        quenchingPerType['redshifts'+str(type)].append(redshifts2_all[j])
        quenchingPerType['mass'+str(type)].append(m)
        quenchingPerType['f_gas'+str(type)].append(frac_gas2_all[j])
    edgcolors = ['b', 'r', 'g']
    mass_ranges = [9.5,10.3,11.0,18.0]
    types = [r'$9.5\leq \log(M_*) < 10.3$', r'$10.3\leq \log(M_*) < 11.0$', r'$\log(M_*) \geq 11.0$']
    for i in range(0, len(mass_ranges)-1):
        red_cent,frequency,frequency_sig,times,times_sig = quenching_histogram(results_folder+'/redshifts_m100n1024.txt',galaxies,max_ngal,mass_ranges[i],mass_ranges[i+1],quenchingPerType['times'+str(i)],
                                                                                quenchingPerType['redshifts'+str(i)], 10)
        axes[0].errorbar(red_cent, frequency,yerr=frequency_sig, label=types[i], marker='o', linestyle='--', capsize=3, markersize=8)
        axes[1].errorbar(red_cent, times, yerr = times_sig, label=types[i], marker='o', linestyle='--', capsize=3, markersize=8)
    axes[0].set_ylabel('Quenching events per galaxy', fontsize=16)
    axes[1].set_ylabel(r' $\langle \log(t_{q}/t_{U}) \rangle$', fontsize=16)
    axes[1].set_xlabel('z', fontsize=16)
    axes[0].legend(loc='best', prop={'size': 12})
    fig2.tight_layout()
    fig2.savefig(str(results_folder)+'quenching_histograms.png', format='png', dpi=200)

print('-----------------------------------------------------')
print('QUENCHING AND REJUVENATION PLOT UTILITIES')
print(' ')
print('The following functions are available:')
print(' ')
print('- Scatter plot showing the duration of the different quenching types of events wrt redshift and mass. (Press 1)')
print(' ')
print('- Rejuvenation rate evolution with redshift. (Press 2)')
print(' ')
print('- Quenching times histogram for the frequency and time evolution with redshift for three mass bins. (Press 3)')
print(' ')
u_selec = input('Write the number of the function you would like to use: ')
if u_selec==1:
    Quenching_Scatter_Plot(redshifts2, quenching_times2, ste_mass2)
elif u_selec==2:
    Rejuvenation_Rate_Plot(d, reju_z, counterfile, timefile, redshiftfile)
elif u_selec==3:
    Quenching_Histogram_Plots(quenching_times2_all, redshifts2_all, ste_mass2_all, frac_gas2_all)
else:
    print('ERROR: function not found')
