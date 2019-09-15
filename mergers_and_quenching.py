#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 26 June 2019

This code is an example analysis of the results from the mergerFinder and quenchingFinder code. In this case, the 
analysis performed is the study of the relation between the ocurrence of mergers, quenchings and rejuvenations.

The version here detailed provides the relation plots for mergers and quenching given in Rodriguez et al. (2019).

@author: currorodriguez
"""
# Import required libraries
import numpy as np
import matplotlib
from scipy import stats
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
sns.set(style="white")
import cPickle as pickle
import sys

MODEL = sys.argv[1]  # e.g. m50n512
WIND = sys.argv[2]  # e.g. s50 for Simba

# Import other codes
from quenchingFinder import GalaxyData
results_folder = '../mandq_relations/%s/' % (MODEL) # You can change this to the folder where you want your resulting plots
merger_file = '../mergers/%s/merger_results.pkl' % (MODEL) # File holding the progen info of galaxies
quench_file = '../quench_analysis/%s/quenching_results.pkl' % (MODEL) # File holding the progen info of galaxies

# Extract data from mergers and quenching pickle files
obj = open(merger_file, 'rb')
merger_data = pickle.load(obj)
obj.close()
mergers, sf_galaxies, max_redshift_mergers = merger_data['mergers'], merger_data['sf_galaxies'], merger_data['redshift_limit']

obj = open(quench_file, 'rb')
quench_data = pickle.load(obj)
obj.close()
galaxies_interpolated = quench_data['quenched_galaxies']
mass_limit = quench_data['mass_limit']
# # Save results of rejuvenations coming from first loop
# reju_z = []
# reju_m = []
# reju_t = []
# reju_id = []

# for i in range(len(galaxies_interpolated)):
#     galaxy = galaxies_interpolated[i]
#     for k in range(0, len(galaxy.rate), 4):
#         #if np.log10(galaxy.rate[k+1])>=mass_limit:
#         #print(galaxy.rate)
#         reju_z.append(galaxy.rate[k])
#         reju_t.append(galaxy.rate[k+1])
#         reju_m.append(galaxy.rate[k+2])
#         reju_id.append(galaxy.rate[k+3])

# # Save results of quenchings in the classification scheme chosen
# redshifts2_all = []
# quenching_times2_all = []
# ste_mass2_all = []
# frac_gas2_all = []
# thubble2_all = []

# redshifts2 = [[[],[],[]],[[],[],[]]]
# quenching_times2 = [[[],[],[]],[[],[],[]]]
# ste_mass2 = [[[],[],[]],[[],[],[]]]
# frac_gas2 = [[[],[],[]],[[],[],[]]]
# thubble2 = [[[],[],[]],[[],[],[]]]
# sfr_2 = [[[],[],[]],[[],[],[]]]


# finalis = 0
# nofinalis = 0
# for i in range(0, len(galaxies_interpolated)):
#     galaxy = galaxies_interpolated[i]
#     if len(galaxy.quenching)>0:
#         lastquench = galaxy.quenching[-1]
#     for quench in galaxy.quenching:
#         start = quench.above9 + 1
#         end = quench.below11
#         if quench is lastquench:
#             pos = 0
#             finalis = finalis + 1
#             for k in range(0, len(galaxy.rate), 3):
#                 if 0.2*galaxy.galaxy_t[start]> (galaxy.galaxy_t[start] - galaxy.rate[k+1]) >=0:
#                     pos = 2
#         else:
#             nofinalis = nofinalis + 1
#             pos = 1
#         if np.log10(galaxy.m_gal)>=mass_limit:
#             redshifts2[int(quench.type)][pos].append(galaxy.z_gal)
#             ste_mass2[int(quench.type)][pos].append(np.log10(galaxy.m_gal))
#             quenching_times2[int(quench.type)][pos].append(np.log10(quench.quench_time))
#             frac_gas2[int(quench.type)][pos].append(np.log10(galaxy.fgas_gal))
#             thubble2[int(quench.type)][pos].append(np.log10(galaxy.galaxy_t[end]))
#             sfr_2[int(quench.type)][pos].append(np.log10(galaxy.ssfr_gal[start-1]*galaxy.m_gal))
#             redshifts2_all.append(galaxy.z_gal)
#             ste_mass2_all.append(galaxy.m_gal)
#             quenching_times2_all.append(quench.quench_time)
#             frac_gas2_all.append(galaxy.fgas_gal)
#             thubble2_all.append(galaxy.galaxy_t[end])

def mqr_relation():
    print('Start finding for connection between mergers and quenching')
    time_diff_s = []
    time_diff_f = []
    merger_ratios = []
    quenching_times = []
    for i in range(0, len(mergers)):
        merg = mergers[i]
        possible_q = []
        possible_q_time = []
        possible_r = []
        for j in range(0, len(galaxies_interpolated)):
            galaxy = galaxies_interpolated[j]
            if galaxy.id == merg.id:
                for quench in galaxy.quenching:
                    start = quench.above9 + 1
                    end = quench.below11
                    if np.log10(galaxy.m_gal)>=mass_limit:
                        possible_q.append(galaxy.galaxy_t[start])
                        possible_q_time.append(quench.quench_time/galaxy.galaxy_t[end])
        diff = []
        for k in range(0, len(possible_q)):
            diff.append(possible_q[k] - merg.galaxy_t[1])
        diff = np.asarray(diff)
        if len(possible_q)>0:
            if diff[np.argmin(diff)]>=0 and merg.merger_ratio<=0.6:
                merger_ratios.append(merg.merger_ratio)
                quenching_times.append(possible_q_time[np.argmin(diff)])
                if np.log10(possible_q_time[np.argmin(diff)])>=-1.5:
                    time_diff_s.append(possible_q[np.argmin(diff)]-merg.galaxy_t[1])
                else:
                    time_diff_f.append(possible_q[np.argmin(diff)]-merg.galaxy_t[1])
    time_diff_s = np.asarray(time_diff_s)
    time_diff_f = np.asarray(time_diff_f)
    merger_ratios = np.asarray(merger_ratios)
    quenching_times = np.asarray(quenching_times)
    fig = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r'$t - t_m$(Gyr)', fontsize=16)
    ax.set_ylabel(r'$N/N_{SF}(Total)$', fontsize=16)
    binis = np.linspace(np.log10(1e-1), np.log10(7), 10)
    binis = 10**binis
    time_diff = np.concatenate((time_diff_s,time_diff_f))
    hist, bin_edges = np.histogram(time_diff, bins=binis)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    hist = hist/np.sum(d['sf_galaxies_per_snap'])
    ax.step(bin_cent, hist, 'k', label='Total quenchings', where='mid')
    median = np.median(time_diff)
    ax.plot([median, median],[0, hist.max()], 'k:', linewidth=2.0)
    hist, bin_edges = np.histogram(time_diff_s, bins=binis)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    hist = hist/np.sum(d['sf_galaxies_per_snap'])
    ax.step(bin_cent, hist, 'b', label='Slow quenchings', where='mid')
    median = np.median(time_diff_s)
    ax.plot([median, median],[0, hist.max()], 'b:', linewidth=2.0)
    hist, bin_edges = np.histogram(time_diff_f, bins=binis)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    hist = hist/np.sum(d['sf_galaxies_per_snap'])
    ax.step(bin_cent, hist, 'r', label='Fast quenchings', where='mid')
    median = np.median(time_diff_f)
    ax.plot([median, median],[0, hist.max()], 'r:', linewidth=2.0)


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
                time_diff.append(possible_r[np.argmin(diff)]-merg.galaxy_t[1]+1e-7)
    time_diff = np.asarray(time_diff)
    median = np.median(time_diff)
    merger_boost = np.asarray(merger_boost)
    binis = np.linspace(np.log10(1e-1), np.log10(7), 10)
    binis = 10**binis
    hist, bin_edges = np.histogram(time_diff, bins=binis)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    hist = hist/np.sum(d['sf_galaxies_per_snap'])
    ax.plot([median, median],[0, hist.max()], 'm:', linewidth=2.0)
    ax.step(bin_cent, hist, 'm', label='Rejuvenations', where='mid')
    ax.legend(loc='best', fontsize=16)
    ax.tick_params(labelsize=12)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.set_xlim([0.0,6.66])
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mergertime_and_quench_reju.png',format='png', dpi=250, bbox_inches='tight')

def reju_fastquench(galaxies_interpolated):

    delays = []

    for gal in galaxies_interpolated:
        possible_r = []
        possible_q = []
        for k in range(0, len(gal.rate), 4):
            possible_r.append(gal.rate[k+1])
        for quench in gal.quenching:
            start = quench.above9 + 1
            end = quench.below11
            if np.log10(gal.m_gal[start])>=9.5 and np.log10(quench.quench_time/gal.galaxy_t[end])<-1.5:
                possible_q.append(gal.galaxy_t[start])
        possible_q = np.asarray(possible_q)
        possible_r = np.asarray(possible_r)
        if len(possible_r) and len(possible_q):
            for rej in possible_r:
                diff = possible_q - rej
                if diff[np.argmin(diff)]>=0:
                    delays.append(diff[np.argmin(diff)])
    delays = np.asarray(delays)
    median = np.median(delays)
    print(median)
    binis = np.linspace(np.log10(1e-1), np.log10(7), 10)
    binis = 10**binis
    hist, bin_edges = np.histogram(delays, bins=binis)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    fig = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.plot([median, median],[0, hist.max()], 'm:', linewidth=2.0)
    ax.step(bin_cent, hist, 'm', where='mid')
    ax.legend(loc='best', fontsize=16)
    ax.tick_params(labelsize=12)
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.set_xlim([0.0,6.66])
    fig.tight_layout()
    fig.savefig('rej_fastq_delay.png',format='png', dpi=250, bbox_inches='tight')


def quench_merger_scatter():
    print('Start finding for connection between mergers, quenching and rejuvenations')
    quench_t = []
    merger_t = []
    quench_scale = []
    for i in range(0, len(mergers)):
        merg = mergers[i]
        possible_q = []
        possible_q_t = []
        for j in range(0, len(galaxies_interpolated)):
            galaxy = galaxies_interpolated[j]
            if galaxy.id == merg.id:
                for quench in galaxy.quenching:
                    start = quench.above9 #+ 1
                    end = quench.below11
                    if np.log10(galaxy.m_gal[start])>=mass_limit:
                        possible_q.append(galaxy.galaxy_t[start])
                        possible_q_t.append(np.log10(quench.quench_time/galaxy.galaxy_t[end]))
        diff_q = []
        for k in range(0, len(possible_q)):
            diff_q.append(possible_q[k] - merg.galaxy_t[1])
        diff_q = np.asarray(diff_q)
        if len(possible_q)>0:
            if diff_q[np.argmin(diff_q)]>=0 and merg.merger_ratio<=0.6:
                quench_t.append(possible_q[np.argmin(diff_q)])
                quench_scale.append(possible_q_t[np.argmin(diff_q)])
                merger_t.append(merg.galaxy_t[1])
    quench_t = np.asarray(quench_t)
    quench_scale = np.asarray(quench_scale)
    merger_t = np.asarray(merger_t)
    print('Galaxies with mergers and quenching: '+str(len(merger_t)))
    fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    sc = ax.scatter(quench_t, merger_t, c = quench_scale, cmap='jet', s=15)
    ax.plot([1.9, 14], [1.9, 14], 'k--')
    ax.set_xlim([quench_t.min()*0.9, quench_t.max()*1.1])
    ax.set_ylim([merger_t.min()*0.9, merger_t.max()*1.1])
    ax.set_xlabel(r'$t_q $(Gyr)', fontsize=16)
    ax.set_ylabel(r'$t_m $(Gyr)', fontsize=16)
    cb = fig.colorbar(sc, ax=ax, orientation='horizontal')
    cb.set_label(label=r'$\log(\tau_q/t_H)$', fontsize=16)
    ax.tick_params(labelsize=12)
    fig.tight_layout()
    fig.savefig(str(results_folder)+'merger_quench_scatter.png',format='png', dpi=250, bbox_inches='tight')

def quench_delay(delay,quenching_times, merger_ratios):
    print('Making scatter plot of quenching times vs delay of the quenching process after a merger...')
    fig1 = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig1.add_subplot(1,1,1)
    pear = stats.pearsonr(delay, np.log10(quenching_times))
    print('The Pearson correlation R is given by: '+str(pear))
    sc1 = ax1.scatter(delay, np.log10(quenching_times), c=merger_ratios , cmap='viridis', s=40)
    ax1.set_xlabel(r'$T_q - T_m$', fontsize=16)
    ax1.set_ylabel(r'$\log(t_{quench}/t_{Hubble})$', fontsize=16)
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
            diff.append(possible_r[k] - merg.galaxy_t[2])
        diff = np.asarray(diff)
        if len(possible_r)>0:
            if diff[np.argmin(diff)]>=0:
                if merg.fgas_boost<0:
                    merger_boost.append(0.001)
                else:
                    merger_boost.append(merg.fgas_boost)
                time_diff.append(possible_r[np.argmin(diff)]-merg.galaxy_t[1]+1e-7)
    time_diff = np.asarray(time_diff)
    merger_boost = np.asarray(merger_boost)
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(r'$T_r - T_m$(Gyr)', fontsize=16)
    ax.set_ylabel(r'$\log(N) $(Gyr)', fontsize=16)
    ax.hist(time_diff, bins=12, histtype='step', log=True, color='k')
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mergertime_and_rejuvenation.png',format='png', dpi=250, bbox_inches='tight')
def merger_reju_scatter():
    print('Start finding for connection between mergers and rejuvenations')
    merger_t = []
    rejuvenation_t = []
    merger_boost = []
    for i in range(0, len(mergers)):
        merg = mergers[i]
        possible_r = []
        for j in range(0, len(reju_t)):
            if reju_id[j]==merg.id:
                possible_r.append(reju_t[j])
        diff = []
        for k in range(0, len(possible_r)):
            diff.append(possible_r[k] - merg.galaxy_t[2])
        diff = np.asarray(diff)
        if len(possible_r)>0:
            rejuvenation_t.append(possible_r[np.argmin(diff)])
            merger_t.append(merg.galaxy_t[2])
            b = merg.ssfr_gal[2]/merg.ssfr_gal[1]
            merger_boost.append(b)
    merger_t = np.asarray(merger_t)
    rejuvenation_t = np.asarray(rejuvenation_t)
    merger_boost = np.asarray(merger_boost)
    fig = plt.figure(num=None, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    pear = stats.pearsonr(rejuvenation_t, merger_t)
    print('The Pearson correlation R is given by: '+str(pear))
    ax.hexbin(rejuvenation_t, merger_t, gridsize=20, cmap='BuGn')
    sc = ax.scatter(rejuvenation_t, merger_t, c = np.log10(merger_boost), cmap='winter')
    ax.plot([3, 13], [3, 13], 'k--')
    ax.set_xlim([rejuvenation_t.min(), rejuvenation_t.max()])
    ax.set_ylim([merger_t.min(), merger_t.max()])
    ax.set_xlabel(r'$T_r $(Gyr)', fontsize=16)
    ax.set_ylabel(r'$T_m $(Gyr)', fontsize=16)
    ax.tick_params(labelsize=12)
    cb = fig.colorbar(sc, ax=ax, orientation='horizontal')
    cb.set_label(label=r'$\log(\Delta sSFR)$', fontsize=16)
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mergertime_and_rejuvenation_scatter.png',format='png', dpi=250, bbox_inches='tight')
#time_diff, q_times, m_ratios = mergerquench_relation()
#quench_delay(time_diff,q_times,m_ratios)
#merger_reju_relation()
#merger_reju_scatter()
#mqr_relation()
quench_merger_scatter()
#reju_fastquench(galaxies_interpolated)
