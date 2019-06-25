#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import random
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns
sns.set(style="white")
from import_progen import importApp
from quenchingFinder import GalaxyData
from mergerFinder import merger_finder, myrunningmedian
import sys
simfolder = '../progen_analysis/m100n1024'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
simname = 'm100n1024'#input('SIMBA simulation version: ')
results_folder = '../mergers/'+str(simname)+'/'

d, ngal = importApp(str(simfolder))
galaxies = []
for i in range(ngal):
    sfr_gal = d['sfr_gal' + str(i)][::-1]
    sfe_gal = d['sfe_gal' + str(i)][::-1]
    z_gal = d['z_gal' + str(i)][::-1]
    galaxy_t = d['galaxy_t' + str(i)][::-1]
    galaxy_m = d['m_gal'+str(i)][::-1]
    fgas_gal = d['h2_gal'+str(i)][::-1]+0.000001
    gal_type = d['gal_type'+str(i)][::-1]
    galaxy = GalaxyData(i, sfr_gal, sfe_gal, z_gal, galaxy_t, galaxy_m, fgas_gal, gal_type)
    galaxies.append(galaxy)

print('Data extracted from PROGEN files.')
print('Performing search of mergers...')
mergers, sf_galaxies = merger_finder(galaxies, 0.2, 10**9.5, 2.5)
print('Search of mergers finished!!')

def after_before_vs_msqPlots(mergers, sf_galaxies):
    ylabels = [r'$\log(sSFR)$',r'$\log(f_{H_2})$',r'$\log(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    merger_labels = ['Before merger','After merger','MSQ non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    colours = ['b','r','k']
    colour_lines = ['c','m']
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    for i in range(0, len(ylabels)):
        fig, axes = plt.subplots(len(titles), 1, sharex='col', num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
        for m in range(0, len(titles)):
            axes[m].set_ylabel(ylabels[i], fontsize=16)
            a = [[],[]]
            b = [[],[]]
            #a_m = []
            #b_m = []
            for j in range(0, len(mergers)):
                if zlimits[m][0] <= mergers[j].z_gal[1] < zlimits[m][1]:
                    a[0].append(np.log10(mergers[j].m_gal[0]))
                    #a_m.append(np.log10(mergers[j].m_gal[1]))
                    a[1].append(np.log10(mergers[j].m_gal[2]))
                    if i==0:
                        b[0].append(np.log10(mergers[j].ssfr_gal[0]))
                        #b_m.append(np.log10(mergers[j].sfr_gal[1]/mergers[j].m_gal[1]))
                        b[1].append(np.log10(mergers[j].ssfr_gal[2]))
                        for s in range(0, len(b)):
                            if m==0 and b[s][-1]<-11:
                                b[s][-1] = -11 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==1 and b[s][-1]<-10.2:
                                b[s][-1] = -10.2 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==2 and b[s][-1]<-10:
                                b[s][-1] = -10 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                    elif i==1:
                        b[0].append(np.log10(mergers[j].fgas_gal[0]))
                        #b_m.append(np.log10(mergers[j].fgas_gal[1]))
                        b[1].append(np.log10(mergers[j].fgas_gal[2]))
                        for s in range(0, len(b)):
                            if m==0 and b[s][-1]<-2:
                                b[s][-1] = -2 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==1 and b[s][-1]<-1.7:
                                b[s][-1] = -1.7 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==2 and b[s][-1]<-1.5:
                                b[s][-1] = -1.5 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                    elif i==2:
                        b[0].append(np.log10(mergers[j].sfe_gal[0]))
                        #b_m.append(np.log10(mergers[j].sfe_gal[1]))
                        b[1].append(np.log10(mergers[j].sfe_gal[2]))
                        for s in range(0, len(b)):
                            if m==0 and b[s][-1]<-11:
                                b[s][-1] = -11 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==1 and b[s][-1]<-10.5:
                                b[s][-1] = -10.5 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
                            elif m==2 and b[s][-1]<-10.2:
                                b[s][-1] = -10.2 + 0.01*random.randint(0,10)
                                axes[m].arrow(a[s][-1],b[s][-1],0,-0.2, head_width=0.018,
                                                                width=0.005, head_length=0.1,
                                                                color=colours[s])
            for k in range(0, len(a)):
                x,y,ysig = myrunningmedian(np.asarray(a[k]),np.asarray(b[k]),15)
                axes[m].scatter(np.asarray(a[k]),np.asarray(b[k]), color=colours[k], label=merger_labels[k],
                                    marker='.', s=10.0, alpha=0.7)
                axes[m].plot(x, y, color = colour_lines[k], linewidth=2.5)
            a = []
            b = []
            for n in range(0, len(sf_galaxies)):
                if zlimits[m][0] <= sf_galaxies[n].z_gal < zlimits[m][1]:
                    a.append(np.log10(sf_galaxies[n].m_gal))
                    if i==0:
                        b.append(np.log10(sf_galaxies[n].ssfr_gal))
                    elif i==1:
                        b.append(np.log10(sf_galaxies[n].fgas_gal))
                    elif i==2:
                        b.append(np.log10(sf_galaxies[n].sfe_gal))
            x,y,ysig = myrunningmedian(np.asarray(a),np.asarray(b),20)
            axes[m].plot(x, y, color = colours[2], label=merger_labels[2])
            axes[m].fill_between(x, y-ysig, y+ysig, facecolor=colours[2], alpha=0.25)
            axes[m].text(0.05, 0.05, titles[m], transform=axes[m].transAxes, fontsize=14,
                            verticalalignment='bottom', bbox=props)
            axes[m].margins(.2)
            axes[m].set_xlim([9.3,11.9])

        axes[len(titles)-1].set_xlabel(r'$\log(M_{*})$', fontsize=16)

        axes[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0., prop={'size': 13})
        fig.tight_layout()
        fig.savefig(str(results_folder)+'merger_'+str(names[i])+'.png', format='png', dpi=200)

def statsMergers(mergers, sf_galaxies, nbins, printresults = True, plot=False):
    ylabels = [r'$\log(sSFR)$',r'$\log(F_{H_2})$',r'$\log(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    stats_results = {}
    for m in range(0, len(names)):
        stats_results[names[m]] = {}
        for i in range(0, len(titles)):
            stats_results[names[m]][titles[i]] = {}
            aft_m = []
            aft = []
            bef_m = []
            bef = []
            msq_m = []
            msq = []
            for j in range(0, len(mergers)):
                if zlimits[i][0] <= mergers[j].z_gal[1] < zlimits[i][1]:
                    bef_m.append(np.log10(mergers[j].m_gal[0]))
                    aft_m.append(np.log10(mergers[j].m_gal[2]))
                    if m==0:
                        bef.append(np.log10(mergers[j].ssfr_gal[0]))
                        aft.append(np.log10(mergers[j].ssfr_gal[2]))
                    elif m==1:
                        bef.append(np.log10(mergers[j].fgas_gal[0]))
                        aft.append(np.log10(mergers[j].fgas_gal[2]))
                    elif m==2:
                        bef.append(np.log10(mergers[j].sfe_gal[0]))
                        aft.append(np.log10(mergers[j].sfe_gal[2]))
            for n in range(0, len(sf_galaxies)):
                if zlimits[i][0] <= sf_galaxies[n].z_gal < zlimits[i][1]:
                    msq_m.append(np.log10(sf_galaxies[n].m_gal))
                    if m==0:
                        msq.append(np.log10(sf_galaxies[n].ssfr_gal))
                    elif m==1:
                        msq.append(np.log10(sf_galaxies[n].fgas_gal))
                    elif m==2:
                        msq.append(np.log10(sf_galaxies[n].sfe_gal))
            aft_m = np.asarray(aft_m)
            aft = np.asarray(aft)
            bef_m = np.asarray(bef_m)
            bef = np.asarray(bef)
            msq_m = np.asarray(msq_m)
            msq = np.asarray(msq)
            maxs = np.array([aft_m.max(),bef_m.max()])
            mins = np.array([aft_m.min(),bef_m.min()])
            bins = np.linspace(mins.min(), maxs.max(), nbins)
            delta = bins[1] - bins[0]
            bin_cent = bins - delta/2
            stats_results[names[m]][titles[i]]['merger_pvalue'] = np.zeros(len(bins)-1)
            stats_results[names[m]][titles[i]]['aftvsbef_pvalue'] = np.zeros(len(bins)-1)
            stats_results[names[m]][titles[i]]['bin_cent'] = np.delete(bin_cent, 0)
            digi_aft = np.digitize(aft_m, bins, right=True)
            digi_bef = np.digitize(bef_m, bins, right=True)
            digi_msq = np.digitize(msq_m, bins, right=True)
            for j in range(1, len(bins)):
                aft_bin = []
                bef_bin = []
                msq_bin = []
                for k in range(0, len(aft_m)):
                    if digi_aft[k]==j:
                        aft_bin.append(aft[k])
                for k in range(0, len(bef_m)):
                    if digi_bef[k]==j:
                        bef_bin.append(bef[k])
                for k in range(0, len(msq_m)):
                    if digi_msq[k]==j:
                        msq_bin.append(msq[k])
                aft_bin = np.asarray(aft_bin)
                bef_bin = np.asarray(bef_bin)
                merger_bin = np.concatenate((aft_bin, bef_bin), axis=None)
                msq_bin = np.asarray(msq_bin)
                statsKS, pvalue = stats.ks_2samp(merger_bin, msq_bin)
                stats_results[names[m]][titles[i]]['merger_pvalue'][j-1] = pvalue
                statsKS, pvalue = stats.ks_2samp(aft_bin, bef_bin)
                stats_results[names[m]][titles[i]]['aftvsbef_pvalue'][j-1] = pvalue
        if printresults==True:
            print('#########################################')
            print('VARIABLE STUDIED: '+str(names[m]))
            print('Statistical significance of merger with respect to star-forming main sequence')
            for w in range(0, len(titles)):
                print('----------------------------------')
                print('Redshift bin considered: '+str(titles[w]))
                for v in range(0, nbins-1):
                    print('................')
                    print('Mass bin center: '+str(stats_results[names[m]][titles[w]]['bin_cent'][v]))
                    print('p-value from KS 2-sample test: '+str(stats_results[names[m]][titles[w]]['merger_pvalue'][v]))
            print('-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.')
            print('Statistical significance of difference between after and before merger data')
            for w in range(0, len(titles)):
                print('----------------------------------')
                print('Redshift bin considered: '+str(titles[w]))
                for v in range(0, nbins-1):
                    print('................')
                    print('Mass bin center: '+str(stats_results[names[m]][titles[w]]['bin_cent'][v]))
                    print('p-value from KS 2-sample test: '+str(stats_results[names[m]][titles[w]]['aftvsbef_pvalue'][v]))

def distanceMSQ(mergers, sf_galaxies, nbins):
    ylabels = [r'$\Delta_{MSQ}(sSFR)$',r'$\Delta_{MSQ}(f_{H_2})$',r'$\Delta_{MSQ}(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    merger_labels = ['After merger','Before merger','Non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    colours = ['r', 'g']
    lines = ['-','--','-.']
    markers = ['o','v', 's']
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    distance_results = {}
    for i in range(0, len(ylabels)):
        distance_results[names[i]] = {}
        fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
        ax = fig.add_subplot(1,1,1)
        ax.set_ylabel(ylabels[i], fontsize=16)
        ax.set_xlabel(r'$\log(M_{*})$', fontsize=16)
        for m in range(0, len(titles)):
            distance_results[names[i]]['aft_d'+str(m)] = []
            distance_results[names[i]]['aft_cent'+str(m)] = []
            distance_results[names[i]]['bef_d'+str(m)] = []
            distance_results[names[i]]['bef_cent'+str(m)] = []
            aft_m = []
            aft = []
            bef_m = []
            bef = []
            msq_m = []
            msq = []
            for j in range(0, len(mergers)):
                if zlimits[i][0] <= mergers[j].z_gal[1] < zlimits[i][1]:
                    bef_m.append(np.log10(mergers[j].m_gal[0]))
                    aft_m.append(np.log10(mergers[j].m_gal[2]))
                    if m==0:
                        bef.append(np.log10(mergers[j].ssfr_gal[0]))
                        aft.append(np.log10(mergers[j].ssfr_gal[2]))
                    elif m==1:
                        bef.append(np.log10(mergers[j].fgas_gal[0]))
                        aft.append(np.log10(mergers[j].fgas_gal[2]))
                    elif m==2:
                        bef.append(np.log10(mergers[j].sfe_gal[0]))
                        aft.append(np.log10(mergers[j].sfe_gal[2]))
            for n in range(0, len(sf_galaxies)):
                if zlimits[i][0] <= sf_galaxies[n].z_gal < zlimits[i][1]:
                    msq_m.append(np.log10(sf_galaxies[n].m_gal))
                    if m==0:
                        msq.append(np.log10(sf_galaxies[n].ssfr_gal))
                    elif m==1:
                        msq.append(np.log10(sf_galaxies[n].fgas_gal))
                    elif m==2:
                        msq.append(np.log10(sf_galaxies[n].sfe_gal))
            aft_m = np.asarray(aft_m)
            aft = np.asarray(aft)
            bef_m = np.asarray(bef_m)
            bef = np.asarray(bef)
            msq_m = np.asarray(msq_m)
            msq = np.asarray(msq)
            maxs = np.array([aft_m.max(),bef_m.max(),msq_m.max()])
            mins = np.array([aft_m.min(),bef_m.min(),msq_m.min()])
            bins = np.linspace(mins.min(), maxs.max(), nbins)
            delta = bins[1] - bins[0]
            bin_cent = bins - delta/2
            idx = np.digitize(aft_m, bins)
            running_median = [np.median(aft[idx==k]) for k in range(0,nbins)]
            aft_median = np.asarray(running_median)
            idx = np.digitize(bef_m, bins)
            running_median = [np.median(bef[idx==k]) for k in range(0,nbins)]
            bef_median = np.asarray(running_median)
            idx = np.digitize(msq_m, bins)
            running_median = [np.median(msq[idx==k]) for k in range(0,nbins)]
            msq_median = np.asarray(running_median)
            for j in range(0, len(aft_median)):
                if not np.isnan(aft_median[j]) and not np.isnan(msq_median[j]):
                    distance_results[names[i]]['aft_d'+str(m)].append((10**aft_median[j]-10**msq_median[j])/abs(10**msq_median[j]))
                    distance_results[names[i]]['aft_cent'+str(m)].append(bin_cent[j])

            for j in range(0, len(bef_median)):
                if not np.isnan(bef_median[j]) and not np.isnan(msq_median[j]):
                    distance_results[names[i]]['bef_d'+str(m)].append((10**bef_median[j]-10**msq_median[j])/abs(10**msq_median[j]))
                    distance_results[names[i]]['bef_cent'+str(m)].append(bin_cent[j])
            # if names[i]=='gas_frac':
            #     prob = ['aft', 'bef']
            #     for l in range(0, 2):
            #         for k in range(0, len(distance_results[names[i]][str(prob[l])+'_cent'+str(m)])):
            #             if abs(distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]) > 2.5:
            #                 if distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]>0:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = 2.5
            #                     head = 2.5*0.1
            #                 else:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = -2.5
            #                     head = -2.5*0.1
            #                 ax.arrow(distance_results[names[i]][str(prob[l])+'_cent'+str(m)][k],distance_results[names[i]][str(prob[l])+'_d'+str(m)][k],0,head,
            #                             head_width=0.015, width=0.008,
            #                             head_length=0.1, color=colours[l])
            # elif names[i]=='sfe_gal':
            #     prob = ['aft', 'bef']
            #     for l in range(0, 2):
            #         for k in range(0, len(distance_results[names[i]][str(prob[l])+'_cent'+str(m)])):
            #             if abs(distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]) > 0.1:
            #                 if distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]>0:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = 0.1
            #                     head = 0.1*0.1
            #                 else:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = -0.1
            #                     head = -0.1*0.1
            #                 ax.arrow(distance_results[names[i]][str(prob[l])+'_cent'+str(m)][k],distance_results[names[i]][str(prob[l])+'_d'+str(m)][k],0,head,
            #                             head_width=0.015, width=0.008,
            #                             head_length=0.005, color=colours[l])
            # elif names[i]=='burst_ssfr':
            #     prob = ['aft', 'bef']
            #     for l in range(0, 2):
            #         for k in range(0, len(distance_results[names[i]][str(prob[l])+'_cent'+str(m)])):
            #             if abs(distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]) > 0.1:
            #                 if distance_results[names[i]][str(prob[l])+'_d'+str(m)][k]>0:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = 0.1
            #                     head = 0.1*0.1
            #                 else:
            #                     distance_results[names[i]][str(prob[l])+'_d'+str(m)][k] = -0.1
            #                     head = -0.1*0.1
            #                 ax.arrow(distance_results[names[i]][str(prob[l])+'_cent'+str(m)][k],distance_results[names[i]][str(prob[l])+'_d'+str(m)][k],0,head,
            #                             head_width=0.015, width=0.008,
            #                             head_length=0.01, color=colours[l])

            ax.plot(distance_results[names[i]]['aft_cent'+str(m)], distance_results[names[i]]['aft_d'+str(m)], label='After at '+str(titles[m]), color=colours[0], linestyle=lines[m], marker=markers[m])
            ax.plot(distance_results[names[i]]['bef_cent'+str(m)], distance_results[names[i]]['bef_d'+str(m)], label='Before at '+str(titles[m]), color=colours[1], linestyle=lines[m], marker=markers[m])
        ax.plot([9.5, 11.8],[0.0, 0.0], 'k--')
        ax.legend(loc='best')
        #ax.margins(x=None, y=.2)
        fig.tight_layout()
        fig.savefig(str(results_folder)+'distance_msq_'+str(names[i])+'.png', format='png', dpi=200)

def distanceMSQ_2(mergers, sf_galaxies, nbins):
    ylabels = [r'$\Delta_{MSQ}(sSFR)$',r'$\Delta_{MSQ}(f_{H_2})$',r'$\Delta_{MSQ}(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    colours = ['b','r','m']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    fig1 = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
    ax1 = fig1.add_subplot(1,1,1)
    fig2, axes = plt.subplots(3, 1, sharex=True, num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
    fig2.subplots_adjust(hspace=0)
    axes[2].set_xlabel(r'$\log(M_{*})$', fontsize=16)
    bins = np.linspace(9.5, 12.0, nbins)
    delta = bins[1] - bins[0]
    bin_cent = bins - delta/2
    fgas = [0,0,0]
    sfe = [0,0,0]
    for i in range(0, len(ylabels)):
        axes[i].set_ylabel(ylabels[i], fontsize=16)
        axes[i].plot([9.5,12.0],[0.0,0.0], 'k--')
        for m in range(0, len(titles)):
            mer_m = []
            mer = []
            msq_m = []
            msq = []
            for j in range(0, len(mergers)):
                if zlimits[m][0] <= mergers[j].z_gal[1] < zlimits[m][1]:
                    mer_m.append(np.log10(mergers[j].m_gal[1]))
                    if i==0:
                        mer.append(mergers[j].ssfr_gal[1])
                    elif i==1:
                        mer.append(mergers[j].fgas_gal[1])
                    elif i==2:
                        mer.append(mergers[j].sfe_gal[1])
            for n in range(0, len(sf_galaxies)):
                if zlimits[m][0] <= sf_galaxies[n].z_gal < zlimits[m][1]:
                    msq_m.append(np.log10(sf_galaxies[n].m_gal))
                    if i==0:
                        msq.append(sf_galaxies[n].ssfr_gal)
                    elif i==1:
                        msq.append(sf_galaxies[n].fgas_gal)
                    elif i==2:
                        msq.append(sf_galaxies[n].sfe_gal)
            mer_m = np.asarray(mer_m)
            mer = np.asarray(mer)
            msq_m = np.asarray(msq_m)
            msq = np.asarray(msq)
            idx = np.digitize(mer_m, bins)
            running_median = [np.median(mer[idx==k]) for k in range(0,nbins)]
            running_std_1 = np.asarray([np.std(mer[idx==k])/np.sqrt(len(mer[idx==k])) for k in range(0,nbins)])
            mer_median = np.asarray(running_median)
            idx = np.digitize(msq_m, bins)
            running_median = [np.median(msq[idx==k]) for k in range(0,nbins)]
            running_std_2 = np.asarray([np.std(msq[idx==k])/np.sqrt(len(msq[idx==k])) for k in range(0,nbins)])
            msq_median = np.asarray(running_median)
            print(running_std_1,running_std_2)
            distance_std = np.sqrt(((1/msq_median)**2)*running_std_1**2 + (((mer_median/(msq_median**2)-1)**2)*running_std_2**2))
            distance = (mer_median-msq_median)/msq_median
            if i==1:
                fgas[m] = distance
            elif i==2:
                sfe[m] = distance
            axes[i].plot(bin_cent, distance, label=titles[m], color=colours[m])
            axes[i].fill_between(bin_cent, distance-distance_std, distance+distance_std, facecolor=colours[m], alpha=0.25)
            axes[i].set_ylabel(ylabels[i], fontsize=16)
    axes[0].legend(loc='best', prop={'size': 12})
    fig2.savefig(str(results_folder)+'distance_msq.png', format='png', dpi=200)
    for i in range(0, len(fgas)):
        ax1.plot(bin_cent, fgas[i]*sfe[i], label=titles[i])
    ax1.legend(loc='best', prop={'size': 12})
    fig1.savefig(str(results_folder)+'test_distance.png', format='png', dpi=200)
print('MERGER-INDUCED STARBURST ANALYSIS')
print(' ')
print('---------------------------------')
print(' ')
print('The following functions are available:')
print(' ')
print('- Comparison plots showing the after and before of mergers with respect to MSQ. (Press 1)')
print(' ')
print('- 2-sample Kolmogorov-Smirnov test comparing the difference of the after and before merger values with the MSQ. (Press 2)')
print(' ')
print('- Deviation plot of running median of after and before merger with respect to MSQ. (Press 3)')
print(' ')
print('- Deviation plot of running median with respect to MSQ. (Press 4)')
print(' ')
u_selec = input('Write the number of the function you would like to use: ')
if u_selec==1:
    after_before_vs_msqPlots(mergers, sf_galaxies)
elif u_selec==2:
    statsMergers(mergers, sf_galaxies, 5)
elif u_selec==3:
    distanceMSQ(mergers, sf_galaxies, 10)
elif u_selec==4:
    distanceMSQ_2(mergers, sf_galaxies, 10)
else:
    print('ERROR: function not found')
