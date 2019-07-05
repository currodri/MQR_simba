#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
from decimal import Decimal
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns
sns.set(style="ticks")
from import_progen import importApp
from quenchingFinder import GalaxyData
from mergerFinder import merger_finder, myrunningmedian, plotmedian
import sys
import pickle
simfolder = '../progen_analysis/m100n1024'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../mergers/'+str(simname)+'/'
pickle_file = '../progen_analysis/m100n1024/progen_'+str(simname)+'.pkl'

#d, ngal = importApp(str(simfolder))
obj = open(pickle_file, 'rb')
d = pickle.load(obj)
ngal = 49215
galaxies = []
for i in range(ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['t_gal' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]#+0.00000001
    gal_type = d['gal_type'+str(i)][::-1]
    gal_pos = d['gal_pos'+str(i)][::-1]
    galaxy = GalaxyData(i, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type,gal_pos)
    galaxies.append(galaxy)

mergers, sf_galaxies = merger_finder(galaxies, 0.2, 10**9.5, 3.5)

def SF_Budget(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    f_budget = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        sfr_m = 0
        sfr_nm = 0
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                sfr_m = sfr_m + merger.sfr_gal[1]
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                sfr_nm = sfr_nm + msq.ssfr_gal
        if sfr_m != 0 and sfr_nm != 0:
            f_budget[i] = sfr_m/(sfr_m+sfr_nm)
    plt.plot(z_cent, f_budget, linestyle='--', marker='o')
    plt.xlabel(r'$z$')
    plt.ylabel('Fraction of SF Budget')
    plt.savefig(str(results_folder)+'sfbudget.png', dpi=250)

def SFR_Evolution(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    sfr_m_ave = np.zeros(n_bins-1)
    sfr_m_error = np.zeros(n_bins-1)
    sfr_nm_ave = np.zeros(n_bins-1)
    sfr_nm_error = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        sfr_m = []
        sfr_nm = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                sfr_m.append(np.log10(merger.sfr_gal[1]/merger.m_gal[1]))
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                sfr_nm.append(np.log10(msq.ssfr_gal/msq.m_gal))
        sfr_m = np.asarray(sfr_m)
        sfr_nm = np.asarray(sfr_nm)
        sfr_m_ave[i] = np.average(sfr_m)
        sfr_m_error[i] = np.std(sfr_m)/np.sqrt(np.log10(len(sfr_m)))
        sfr_nm_ave[i] = np.average(sfr_nm)
        sfr_nm_error[i] = np.std(sfr_nm)/np.sqrt(np.log10(len(sfr_nm)))
    plt.errorbar(z_cent, sfr_m_ave, yerr=sfr_m_error, linestyle='--', marker='o', label='Mergers star-forming', capsize=2, capthick=2)
    plt.errorbar(z_cent, sfr_nm_ave, yerr=sfr_nm_error, linestyle='--', marker='s', label='Non-mergers star-forming ', capsize=2, capthick=2)
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\log(\langle$sSFR (yr$^{-1}\rangle)$')
    plt.tight_layout()
    plt.legend(loc='best')
    plt.savefig(str(results_folder)+'sfr_evolution.png', dpi=250)
def SFR_Evolution2(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    ssfr_m_ave = np.zeros(n_bins-1)
    ssfr_m_error = np.zeros(n_bins-1)
    ssfr_nm_ave = np.zeros(n_bins-1)
    ssfr_nm_error = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    ssfr_m = []
    ssfr_nm = []
    pos_nm = []
    pos_m = []
    red_m = []
    red_nm = []
    for i in range(0, n_bins-1):
        mergers_m = []
        msq_m = []
        msq_idx = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                ssfr_m.append(np.log10(merger.sfr_gal[1]/merger.m_gal[1]))
                pos_m.append(merger.gal_pos[1])
                mergers_m.append(np.log10(merger.m_gal[1]))
                red_m.append(merger.z_gal[1])
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                msq_idx.append(k)
                msq_m.append(np.log10(msq.m_gal))
        mergers_m = np.asarray(mergers_m)
        msq_idx = np.asarray(msq_idx)
        msq_m = np.asarray(msq_m)
        msq_sort_idx = np.argsort(msq_m)
        msq_sorted = np.sort(msq_m)
        msq_sel = []
        for m in range(0, len(mergers_m)):
            if len(msq_sorted)>3:
                loc = np.searchsorted(msq_sorted, mergers_m[m])
                below = 30
                above = 30
                if loc - below < 0:
                    below = loc
                if loc + above > len(msq_sorted):
                    above  = len(msq_sorted) - loc
                msq_slice = msq_sorted[loc-below:loc+above]
                msq_diff = msq_slice - mergers_m[m]
                for difmin in range(0, 8):
                    min_idx = np.argmin(msq_diff)
                    idx = loc - below + min_idx
                    msq_sel.append(msq_sort_idx[idx])
                    msq_diff = np.delete(msq_diff, min_idx)
                    msq_sorted = np.delete(msq_sorted, idx)
                    msq_sort_idx = np.delete(msq_sort_idx, idx)
        for selected  in range(0, len(msq_sel)):
            real_idx = msq_idx[msq_sel[selected]]
            msq_gal = msq_galaxies[real_idx]
            ssfr_nm.append(np.log10(msq_gal.ssfr_gal/msq_gal.m_gal))
            pos_nm.append(msq.gal_pos)
            red_nm.append(msq.z_gal)
    ssfr_m = np.asarray(ssfr_m)
    ssfr_nm = np.asarray(ssfr_nm)
    pos_m = np.asarray(pos_m)
    pos_nm = np.asarray(pos_nm)
    red_m = np.asarray(red_m)
    red_nm = np.asarray(red_nm)
    print(pos_m)
    print(pos_nm)
    print(red_m)
    print(red_nm)

    ssfr_m_ave, ssfr_m_error = plotmedian(red_m,ssfr_m, pos=pos_m, boxsize=d['boxsize_in_kpccm'],bin_choosen=z_bins)
    ssfr_nm_ave, ssfr_nm_error = plotmedian(red_nm,ssfr_nm, pos=pos_nm, boxsize=d['boxsize_in_kpccm'],bin_choosen=z_bins)
    return ssfr_m_ave,ssfr_m_error,ssfr_nm_ave,ssfr_nm_error,z_cent
    # plt.errorbar(z_cent, ssfr_m_ave, yerr=ssfr_m_error, linestyle='--', marker='o', label='Mergers star-forming', capsize=2, capthick=2)
    # plt.errorbar(z_cent, ssfr_nm_ave, yerr=ssfr_nm_error, linestyle='--', marker='s', label='Mass-matched sample of non-merger star-forming', capsize=2, capthick=2)
    # plt.xlabel(r'$z$')
    # plt.ylabel(r'$\log(\langle$sSFR (yr$^{-1}\rangle)$')
    # plt.tight_layout()
    # plt.legend(loc='best')
    # plt.savefig(str(results_folder)+'sfr_evolution_2.eps', dpi=250)
def Merger_Fraction(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    f_merger = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        m_counter = 0
        nm_counter = 0
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                m_counter = m_counter + 1
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                nm_counter = nm_counter + 1
        f_merger[i] = m_counter/(m_counter+nm_counter)
    plt.plot(z_cent, f_merger, linestyle='--', marker='o')
    plt.xlabel(r'$z$')
    plt.ylabel('Merger fraction of star-forming galaxies')
    plt.savefig(str(results_folder)+'mfr_evolution.png', dpi=250)
def Merger_Fraction_Mass_Distribution(mergers, msq_galaxies, n_bins):
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    markers = ['o','v', 's']
    mass_bins = np.linspace(9.5, 12, n_bins)
    delta = mass_bins[1]-mass_bins[0]
    mass_cent = mass_bins - delta/2
    mass_cent = np.delete(mass_cent, 0)
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    for zs in range(0, len(zlimits)):
        f_merger = np.zeros(n_bins-1)
        for i in range(0, n_bins-1):
            m_counter = 0
            nm_counter = 0
            for j in range(0, len(mergers)):
                merger = mergers[j]
                if zlimits[zs][0] <= merger.z_gal[1] < zlimits[zs][1] and mass_bins[i] <= np.log10(merger.m_gal[1]) < mass_bins[i+1]:
                    m_counter = m_counter + 1
            for k in range(0, len(msq_galaxies)):
                msq = msq_galaxies[k]
                if zlimits[zs][0] <= msq.z_gal < zlimits[zs][1] and mass_bins[i] <= np.log10(msq.m_gal) < mass_bins[i+1]:
                    nm_counter = nm_counter + 1
            if m_counter != 0 or nm_counter != 0:
                f_merger[i] = float(m_counter)/(float(m_counter)+float(nm_counter))
        ax.plot(mass_cent, f_merger, label=titles[zs], linestyle='--', marker=markers[zs])
    ax.set_xlabel(r'$\log(M_{*})$', fontsize=16)
    ax.legend(loc='best')
    ax.set_ylabel('Merger fraction of star-forming galaxies')
    fig.tight_layout()
    fig.savefig(str(results_folder)+'mfr_evolution_permass.png', dpi=250)
def Merger_Contribution(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    f_merger = np.zeros(n_bins-1)
    f_budget = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        m_counter = 0
        nm_counter = 0
        sfr_m = 0
        sfr_nm = 0
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                sfr_m = sfr_m + merger.sfr_gal[1]
                m_counter = m_counter + 1
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                sfr_nm = sfr_nm + msq.ssfr_gal
                nm_counter = nm_counter + 1
        f_merger[i] = float(Decimal(m_counter)/Decimal(m_counter+nm_counter))
        f_budget[i] = sfr_m/(sfr_m+sfr_nm)
    f_merger = np.asarray(f_merger)
    f_budget = np.asarray(f_budget)
    z_cent = np.asarray(z_cent)
    return f_merger,f_budget,z_cent
    # fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    # ax = fig.add_subplot(1,1,1)
    # ax.set_xscale("log")
    # ax.set_xticks([0.0, 0.5, 1.0, 2.5])
    # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    # ax.set_yscale("log")
    # plt.plot(z_cent, f_merger, linestyle='--', marker='o', label='Fraction of galaxies')
    # plt.plot(z_cent, f_budget, linestyle='--', marker='s', label='Fraction of SF Budget')
    # ax.set_xlabel(r'$z$', fontsize=16)
    # ax.set_ylabel('Merger contribution to star-forming galaxies', fontsize=16)
    # ax.legend(loc='best')
    # fig.tight_layout()
    # fig.savefig(str(results_folder)+'merger_contribution.png',format='png', dpi=250)

def Fgas_mean(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    fgas_ave = np.zeros(n_bins-1)
    m_ave = np.zeros(n_bins-1)
    fgas_error = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        fgas = []
        mass = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                fgas.append(merger.fgas_gal[1])
                mass.append(merger.m_gal[1])
        mass = np.asarray(mass)
        m_ave[i] = np.average(mass)
        fgas = np.asarray(fgas)
        fgas_ave[i] = np.average(fgas)
        fgas_error[i] = np.std(fgas)/np.sqrt(len(fgas))
    model_1 = 0.04*(m_ave/(4.5e+11))**(-0.59*(1+z_cent)**(0.45))
    fig = plt.figure(num=None, figsize=(8, 4), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    ax.set_xscale("log")
    ax.set_xticks([0.0, 0.5, 1.0, 2.5])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yscale("log")
    plt.errorbar(z_cent, fgas_ave, yerr=fgas_error, linestyle='--', marker='o', capsize=2, capthick=2, label=r'Mergers in m50n512')
    plt.plot(z_cent, model_1,  'k-', label='Stewart et al. 2009b')
    plt.errorbar([1.465, 1.523, 1.522, 1.414, 1.6, 1.459],
                    [0.66, 0.51, 0.58, 0.62, 0.52,0.56],
                    yerr=[0.09*0.66, 0.09*0.51, 0.09*0.58, 0.09*0.62, 0.09*0.52, 0.09*0.56],
                    label='Daddi et al. 2010', marker='d', capsize=2, linestyle=' ',markerfacecolor='None')
    plt.errorbar([1.2,1.17,1.18,1.12,1.19,1.23,1.28,1.12,1.1,2.19,2.33,2.43,2.29,2.34,2.34,2.19,2.17,2.18,2.21],
                    [0.46,0.34,0.47,0.27,0.45,0.46,0.18,0.14,0.26,0.51,0.73,0.42,0.36,0.78,0.4,0.32,0.33,0.41,0.57],
                    yerr=[0.19,0.14,0.19,0.11,0.2,0.19,0.08,0.07,0.12,0.26,0.32,0.22,0.15,0.33,0.17,0.16,0.17,0.18,0.23],
                    label='Tacconi et al. 2010', marker='s', capsize=2, linestyle=' ',markerfacecolor='None')
    ax.set_xlabel(r'$z$', fontsize=16)
    ax.set_ylabel(r'$\langle f_{H_2}\rangle$', fontsize=16)
    ax.legend(loc='best')
    fig.tight_layout()
    fig.savefig(str(results_folder)+'fgas_evolution.png',format='png', dpi=250)

def Frac_Merger_rate(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    f_merger = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        m_counter = 0
        nm_counter = 0
        times = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                m_counter = m_counter + 1
                times.append(merger.galaxy_t[1])
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                nm_counter = nm_counter + 1
        times = np.asarray(times)
        delta_t = times.max() - times.min()
        f_merger[i] = m_counter/((m_counter+nm_counter)*delta_t)
    z_cent = np.asarray(z_cent)
    f_merger = np.asarray(f_merger)
    x = np.log10(1+z_cent)
    y = np.log10(f_merger)
    def func(x, a, b):
        return a*x + b
    # popt, pcov = curve_fit(func, x, y)
    # perr = np.sqrt(np.diag(pcov))
    # print(popt)
    # print(perr)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("slope: %f    intercept: %f    r_value: %f    p_value: %f    std_error: %f" % (slope, intercept,r_value, p_value, std_err))
    fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1,1,1)
    #ax.set_xscale("log")
    #ax.set_xticks([0.0, 0.5, 1.0, 2.5])
    #ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.set_yscale("log")
    plt.plot(np.log10(1+z_cent), np.log10(f_merger), linestyle='--', marker='d', label=r'Mergers in the $100h^{-1}$ Mpc box')
    plt.plot(np.log10(1+z_cent),np.log10((10**intercept)*(1+z_cent)**(slope)), 'k-', label='The best fit power law to this data')
    ax.set_xlabel(r'$\log(1+z)$', fontsize=16)
    ax.set_ylabel(r'$log(\mathcal{R}_{merg})$', fontsize=16)
    ax.legend(loc='best', prop={'size': 12})
    fig.tight_layout()
    fig.savefig(str(results_folder)+'merger_rate_evolution.png',format='png', dpi=250)
def Contribution_and_Rate(mergers, msq_galaxies, n_bins):
    z_bins = np.linspace(0.0, 3.5, n_bins)
    f_merger_n = np.zeros(n_bins-1)
    f_merger = np.zeros(n_bins-1)
    f_budget = np.zeros(n_bins-1)
    delta = z_bins[1]-z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    for i in range(0, n_bins-1):
        m_counter = 0
        nm_counter = 0
        sfr_m = 0
        sfr_nm = 0
        times = []
        for j in range(0, len(mergers)):
            merger = mergers[j]
            if z_bins[i]<= merger.z_gal[1] < z_bins[i+1]:
                m_counter = m_counter + 1
                sfr_m = sfr_m + merger.sfr_gal[1]
                times.append(merger.galaxy_t[1])
        for k in range(0, len(msq_galaxies)):
            msq = msq_galaxies[k]
            if z_bins[i]<= msq.z_gal < z_bins[i+1]:
                sfr_nm = sfr_nm + msq.ssfr_gal
                nm_counter = nm_counter + 1
        f_merger[i] = float(Decimal(m_counter)/Decimal(m_counter+nm_counter))
        f_budget[i] = sfr_m/(sfr_m+sfr_nm)
        times = np.asarray(times)
        delta_t = times.max() - times.min()
        f_merger_n[i] = m_counter/((m_counter+nm_counter)*delta_t)
    z_cent = np.asarray(z_cent)
    f_merger_n = np.asarray(f_merger_n)
    x = np.log10(1+z_cent)
    y = np.log10(f_merger_n)
    def func(x, a, b):
        return a*x + b
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print("slope: %f    intercept: %f    r_value: %f    p_value: %f    std_error: %f" % (slope, intercept,r_value, p_value, std_err))
    fig, ax = plt.subplots(2, 1, sharex='col', num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    ax[0].plot(np.log10(1+z_cent), np.log10(f_merger_n), linestyle='--', marker='d', label=r'Mergers in the $100h^{-1}$ Mpc box')
    ax[0].plot(np.log10(1+z_cent),np.log10((10**intercept)*(1+z_cent)**(slope)), 'k-', label=r'$\log(\mathcal{R}_{merg})=%f\cdot(1+z)+%f$' % (slope, intercept))
    ax[0].set_ylabel(r'$log(\Gamma_{Mer})$', fontsize=16)
    ax[0].legend(loc='best', prop={'size': 12})

    ax[1].plot(np.log10(1+z_cent), np.log10(f_merger), linestyle='--', marker='o', label=r'$\log($Fraction of galaxies)')
    ax[1].plot(np.log10(1+z_cent), np.log10(f_budget), linestyle='--', marker='s', label=r'$\log($Fraction of SF Budget)')
    ax[1].set_xlabel(r'$\log(1+z)$', fontsize=16)
    ax[1].set_ylabel('Merger contribution to star-forming galaxies', fontsize=16)
    ax[1].legend(loc='best')
    fig.tight_layout()
    fig.savefig(str(results_folder)+'merger_contribution_and_rate.png',format='png', dpi=250)

def SFR_Evolution_and_Contribution(mergers, msq_galaxies, n_bins):
    ssfr_m_ave,ssfr_m_error,ssfr_nm_ave,ssfr_nm_error,z_cent1 = SFR_Evolution2(mergers,msq_galaxies,n_bins)
    f_merger,f_budget,z_cent2 = Merger_Contribution(mergers,msq_galaxies,n_bins)
    fig, axes = plt.figure(2, 1, sharex='x',num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
    axes[0].errorbar(np.log10(1+z_cent1), ssfr_m_ave, yerr=ssfr_m_error, linestyle='--', marker='o', label='Mergers star-forming', capsize=2, capthick=2)
    axes[0].errorbar(np.log10(1+z_cent1), ssfr_nm_ave, yerr=ssfr_nm_error, linestyle='--', marker='s', label='Mass-matched sample of non-merger star-forming', capsize=2, capthick=2)
    axes[0].set_ylabel(r'$\log(\langle$sSFR (yr$^{-1}\rangle)$')
    axes[0].legend(loc='best')
    axes[1].plot(np.log10(1+z_cent2), np.log10(f_merger), linestyle='--', marker='o', label=r'$\log($Fraction of galaxies)')
    axes[1].plot(np.log10(1+z_cent2), np.log10(f_budget), linestyle='--', marker='s', label=r'$\log($Fraction of SF Budget)')
    axes[1].set_ylabel('Merger contribution to star-forming galaxies', fontsize=16)
    axes[1].legend(loc='best')
    axes[1].set_xlabel(r'$\log(1+z)$')
    fig.tight_layout()
    fig.savefig(str(results_folder)+'sfr_evo_and_contribution.png',format='png', dpi=250)

#SFR_Evolution(mergers, sf_galaxies, 10)
#Merger_Fraction(mergers, sf_galaxies, 10)
#Fgas_mean(mergers, sf_galaxies, 10)
print(' ')
print(' ')
print('MERGER STATISTICAL ANALYSIS')
print(' ')
print('---------------------------------')
print(' ')
print('The following functions are available:')
print(' ')
print('- Mass-matched comparison of sSFR between mergers and star-forming galaxies. (Press 1)')
print(' ')
print('- Evolution of contribution of mergers to the star-forming population. (Press 2)')
print(' ')
print('- Evolution of merger rate with redshift. (Press 3)')
print(' ')
print('- Evolution of contribution and merger rate with redshift. (Press 4)')
print(' ')
print('- Mass distribution of merger fraction in three redshift bins. (Press 5)')
print(' ')
print('- Mass-matched comparison of sSFR between mergers and star-forming galaxies and contribution. (Press 6)')
print(' ')
print('- If you want to do all of them, just Press 7.')
print(' ')
u_selec = input('Write the number of the function you would like to use: ')

if u_selec==1:
    SFR_Evolution2(mergers, sf_galaxies, 10)
elif u_selec==2:
    Merger_Contribution(mergers, sf_galaxies, 10)
elif u_selec==3:
    Frac_Merger_rate(mergers, sf_galaxies, 15)
elif u_selec==4:
    Contribution_and_Rate(mergers, sf_galaxies, 15)
elif u_selec==5:
    Merger_Fraction_Mass_Distribution(mergers, sf_galaxies, 15)
elif u_selec==6:
    SFR_Evolution_and_Contribution(mergers, sf_galaxies, 9)
elif u_selec==7:
    SFR_Evolution2(mergers, sf_galaxies, 10)
    Merger_Contribution(mergers, sf_galaxies, 10)
    Frac_Merger_rate(mergers, sf_galaxies, 15)
    Contribution_and_Rate(mergers, sf_galaxies, 15)
    Merger_Fraction_Mass_Distribution(mergers, sf_galaxies, 15)
    SFR_Evolution_and_Contribution(mergers, sf_galaxies, 15)
else:
    print('ERROR: function not found')
