#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 18:56:07 2018

This code is an example analysis of the results from the mergerFinder and quenchingFinder code. In this case, the 
analysis performed is classification of quenching events and the analysis of their population distribution.

The version here detailed provides the rate plots given in Rodriguez et al. (2019).
@author: currorodriguez
"""

# Import required libraries
import numpy as np
import cPickle as pickle
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style="white")
import sys

MODEL = sys.argv[1]  # e.g. m50n512
WIND = sys.argv[2]  # e.g. s50 for Simba

# Import other codes
from quenchingFinder import GalaxyData
results_folder = '../quench_analysis/%s/' % (MODEL) # You can change this to the folder where you want your resulting plots
#quench_file = '../quench_analysis/%s/quenching_results.pkl' % (MODEL) # File holding the progen info of galaxies
data_file = '/home/curro/quenchingSIMBA/code/SH_Project/mandq_results_%s.pkl' % (MODEL)

# Extract data from quenching pickle file
# obj = open(quench_file, 'rb')
# quench_data = pickle.load(obj)
# obj.close()
# galaxies_interpolated = quench_data['quenched_galaxies']
print('Loading pickle file with data...')
obj = open(data_file, 'rb')
quench_data = pickle.load(obj)
obj.close()
galaxies_interpolated = quench_data['galaxies']
print('Data extracted from pickle file!')

mass_limit = quench_data['mass_limit']

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
print('Number of quenching events in second loop: '
        +str(sum([1 for galaxy in galaxies_interpolated for quench in galaxy.quenching])))
print('Quenching and Rejuvenation analysis done.')
print(' ')

def Fraction_Fast_vs_Slow(x, times, sf_d, bins):
    slow = []
    fast = []
    delta = bins[1] - bins[0]
    cent = bins - delta/2
    cent = np.delete(cent, 0)
    cent = np.asarray(cent)
    fast_bins = np.zeros(len(bins)-1)
    slow_bins = np.zeros(len(bins)-1)
    fast_counter = 0
    slow_counter = 0
    for i in range(0, len(times)):
        if times[i] < -1.5:
            fast.append(x[i])
        else:
            slow.append(x[i])
    for i in range(0, len(bins)-1):
        sf = 0
        f = 0
        s = 0
        for j in range(0, len(sf_d[0])):
            if len(sf_d)>1:
                if bins[i] <= sf_d[1][j] < bins[i+1]:
                    sf = sf + sf_d[0][j]
            else:
                if bins[i] <= sf_d[0][j] < bins[i+1]:
                    sf = sf + 1
        for j in range(0, len(slow)):
            if bins[i] <= slow[j] < bins[i+1]:
                s = s + 1
        for j in range(0, len(fast)):
            if bins[i] <= fast[j] < bins[i+1]:
                f = f + 1
        fast_counter = fast_counter + f
        slow_counter = slow_counter + s
        fast_bins[i] = float(f)/float(sf)
        slow_bins[i] = float(s)/float(sf)
    print(fast_counter, slow_counter)
    return fast_bins,slow_bins,cent

# Plot the results
def Quenching_Scatter_Plot(redshifts, quenching_times, ste_mass):
    y_labels = [r'$\log(\tau_q/t_H)$',r'$\log(\tau_q/t_H)$',r'$\log(N/N_{SF})$']
    frac_labels = ['Centrals with ','Satellites with ']
    text_labels = ['Central galaxies', 'Satellite galaxies']
    x_labels = [r'$\log(1+z)$', r'$\log(M_*[M_{\odot}])$']
    name_file = ['redshift', 'mass']
    colours = ['r','b']
    linestyles = ['-','--']
    sf_x = [quench_data['redshifts'],np.log10(quench_data['sf_galaxies_mass'])]
    x_data = [redshifts, ste_mass]
    props = dict(boxstyle='round', facecolor='white', edgecolor='k', alpha=0.7)
    for i in range(0, len(x_labels)):
        fig, ax = plt.subplots(3, 1, sharex=True, num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
        figR = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
        axR = figR.add_subplot(1,1,1)
        for j in range(0, len(y_labels)):
            ax[j].set_ylabel(y_labels[j], fontsize=16)
            ax[j].tick_params(labelsize=12)
            if j!=2:
                if i==0:
                    a = np.asarray(x_data[i][j][0])
                    a = np.log10(1+a)
                    b =np.asarray(x_data[i][j][2])
                    b = np.log10(1+b)
                elif i==1:
                    a = np.asarray(x_data[i][j][0])
                    b = np.asarray(x_data[i][j][2])
                ax[j].hexbin(a, quenching_times[j][0], bins='log', cmap='Greys', gridsize=30)
                ax[j].scatter(b, quenching_times[j][2], s=8, alpha=0.8, facecolor='g')
                ax[j].text(0.70, 0.90, text_labels[j], transform=ax[j].transAxes, fontsize=14,
                            verticalalignment='top', bbox=props)
                ax[j].plot([a.min(),a.max()],[-1.5,-1.5], 'k--')
            else:
                for k in range(0, 2):
                    x_datas = x_data[i][k][0] + x_data[i][k][2]
                    x_datas = np.asarray(x_datas)
                    x_datas_r = np.asarray(x_data[i][k][2])
                    quenchs = quenching_times[k][0] + quenching_times[k][2]
                    quenchs = np.asarray(quenchs)
                    quenchs_r = np.asarray(quenching_times[k][2])
                    if i==0:
                        sf_data = [quench_data['sf_galaxies_per_snap'],sf_x[i]]
                        pre_bins = np.linspace(0.001,0.7,12)
                        bins = 10**pre_bins - 1
                    else:
                        sf_data = [sf_x[i]]
                        bins = np.linspace(9.5,12.5,12)
                    fast, slow, cent = Fraction_Fast_vs_Slow(x_datas, quenchs, sf_data, bins)
                    fast_r, slow_r, cent_r = Fraction_Fast_vs_Slow(x_datas_r, quenchs_r, sf_data, bins)
                    if i==0:
                        ax[j].plot(np.log10(1+cent), np.log10(fast), label = frac_labels[k]+'fast quenching', color=colours[0], ls=linestyles[k])
                        ax[j].plot(np.log10(1+cent), np.log10(slow), label = frac_labels[k]+'slow quenching', color=colours[1], ls=linestyles[k])
                        axR.plot(np.log10(1+cent), np.log10(fast), color=colours[k], ls='-')
                        axR.plot(np.log10(1+cent), np.log10(slow), color=colours[k], ls='--')
                        axR.plot(np.log10(1+cent_r), np.log10(fast_r), label = frac_labels[k]+'fast quenching with rejuvenation', color=colours[k], ls=':')
                        axR.plot(np.log10(1+cent_r), np.log10(slow_r), label = frac_labels[k]+'slow quenching with rejuvenation', color=colours[k], ls='-.')
                    elif i==1:
                        ax[j].plot(cent, np.log10(fast), label = frac_labels[k]+'fast quenching', color=colours[0], ls=linestyles[k])
                        ax[j].plot(cent, np.log10(slow), label = frac_labels[k]+'slow quenching', color=colours[1], ls=linestyles[k])
                        axR.plot(cent, np.log10(fast), color=colours[k], ls='-')
                        axR.plot(cent, np.log10(slow), color=colours[k], ls='--')
                        axR.plot(cent_r, np.log10(fast_r), label = frac_labels[k]+'fast quenching with rejuvenation', color=colours[k], ls=':')
                        axR.plot(cent_r, np.log10(slow_r), label = frac_labels[k]+'slow quenching with rejuvenation', color=colours[k], ls='-.')

                ax[j].legend(loc='best', prop={'size': 12}, fontsize=14)
        ax[2].set_xlabel(x_labels[i], fontsize=16)
        if i==0:
            axZ = ax[0].twiny()
            maxlz = 0.7
            ax[0].set_xlim(0,maxlz)
            axZ.set_xlim(0,maxlz)
            topticks1 = np.array([0,1,2,3,4])  # desired redshift labels
            topticks2 = np.log10(1+topticks1)  # tick locations in time
            axZ.set_xticklabels(topticks1)
            axZ.set_xticks(topticks2)
            axZ.xaxis.set_ticks_position('top') # set the position of the second x-axis to top
            axZ.xaxis.set_label_position('top') # set the position of the second x-axis to top
            axZ.set_xlabel('z', fontsize=16)
        fig.subplots_adjust(hspace=0)
        fig.savefig(str(results_folder)+'quenching_scatter_'+str(name_file[i])+'.png', format='png', dpi=250, bbox_inches='tight')
        axR.set_xlabel(x_labels[i], fontsize=16)
        axR.set_ylabel(y_labels[2], fontsize=16)
        axR.legend(loc='best', prop={'size': 10})
        figR.savefig(str(results_folder)+'quenching2_rej_'+str(name_file[i])+'.png', format='png', dpi=250, bbox_inches='tight')
def Quenching_Scatter_Plot2(redshifts2, quenching_times2, ste_mass2):
    scatter_labels = [['Final quenching Sat', 'Non-final quenching Sat', 'Final quenching Sat with rejuvenation' ],['Final quenching Central', 'Non-final quenching Central', 'Final quenching Central with rejuvenation']]
    scatter_markers = ['.','*', '.']
    cmaps = ['Reds', 'Blues']
    grids = [30,70]
    #scatter_colours = [['#e90c35','#ec5b75','#e90c35'],['#0c86e4','#8d0ce4','#0c86e4']]
    scatter_colours = [['#e90c35','m','r'],['#0c86e4','c','b']]
    for i in range(0, len(scatter_labels)):
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
    for i in range(len(scatter_labels)-1, -1, -1):
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
