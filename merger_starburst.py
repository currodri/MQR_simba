#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
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
simfolder = '../progen_analysis/m50n512'#input('SIMBA simulation progen folder: ')
sys.path.insert(0, str(simfolder))
simname = 'm50n512'#input('SIMBA simulation version: ')
results_folder = '../mergers/'+str(simname)+'/'

d, ngal = importApp(str(simfolder))
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

mergers, sf_galaxies = merger_finder(galaxies, 0.2, 10**9.5, 2.5)


def after_before_vs_msqPlots(mergers, sf_galaxies):
    ylabels = [r'$\log(sSFR)$',r'$\log(f_{H_2})$',r'$\log(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    merger_labels = ['Before merger','After merger','MSQ non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    colours = ['g','r', 'k']
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    for i in range(0, len(ylabels)):
        fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
        axes = {}
        for m in range(0, len(titles)):
            axes['redbin'+str(m)] = fig.add_subplot(3,1,m+1)
            axes['redbin'+str(m)].set_ylabel(ylabels[i], fontsize=16)
            a = [[],[]]
            b = [[],[]]
            for j in range(0, len(mergers)):
                if zlimits[m][0] <= mergers[j].z_gal[1] < zlimits[m][1]:
                    a[0].append(np.log10(mergers[j].m_gal[0]))
                    a[1].append(np.log10(mergers[j].m_gal[2]))
                    if i==0:
                        b[0].append(np.log10(mergers[j].sfr_gal[0]))
                        b[1].append(np.log10(mergers[j].sfr_gal[2]))
                    elif i==1:
                        b[0].append(np.log10(mergers[j].fgas_gal[0]))
                        b[1].append(np.log10(mergers[j].fgas_gal[2]))
                    elif i==2:
                        b[0].append(np.log10(mergers[j].sfe_gal[0]))
                        b[1].append(np.log10(mergers[j].sfe_gal[2]))
            for k in range(0, len(a)):
                x,y,ysig = myrunningmedian(np.asarray(a[k]),np.asarray(b[k]),15)
                axes['redbin'+str(m)].scatter(np.asarray(a[k]),np.asarray(b[k]), color=colours[k], label=merger_labels[k], marker='.')
                axes['redbin'+str(m)].plot(x, y, color = colours[k], linewidth=2.5)
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
            axes['redbin'+str(m)].plot(x, y, color = colours[2], label=merger_labels[2])
            axes['redbin'+str(m)].fill_between(x, y-ysig, y+ysig, facecolor=colours[2], alpha=0.25)
            axes['redbin'+str(m)].text(0.05, 0.05, titles[m], transform=axes['redbin'+str(m)].transAxes, fontsize=14,
            verticalalignment='bottom', bbox=props)
            axes['redbin'+str(m)].margins(.2)
            axes['redbin'+str(m)].set_xlim([9.3,11.9])

        axes['redbin'+str(len(titles)-1)].set_xlabel(r'$\log(M_{*})$', fontsize=16)

        axes['redbin0'].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
               ncol=3, mode="expand", borderaxespad=0., prop={'size': 13})
        fig.tight_layout()
        fig.savefig(str(results_folder)+'merger_'+str(names[i])+'.png', format='png', dpi=200)

# redshift_bef, redshift_aft, redshift_all, redshift_mer = merger_finder(galaxies, 0.2, 10**9.5)
#
# merger_data = [redshift_aft, redshift_bef, redshift_all]
#
# def after_before_vs_msqPlots(merger_data):
#     ylabels = [r'$\log(sSFR)$',r'$\log(f_{H_2})$',r'$\log(SFE)$']
#     names = ['burst_ssfr','gas_frac','sfe_gal']
#     merger_labels = ['After merger','Before merger','Non merger']
#     titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
#     colours = ['r', 'g', 'k']
#     props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
#     fig2 = plt.figure(num=None, figsize=(8, 5), dpi=80, facecolor='w', edgecolor='k')
#     ax2 = fig2.add_subplot(1,1,1)
#     for i in range(0, len(ylabels)):
#         fig = plt.figure(num=None, figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
#         axes = {}
#         for m in range(0, len(titles)):
#             axes['redbin'+str(m)] = fig.add_subplot(3,1,m+1)
#             axes['redbin'+str(m)].set_ylabel(ylabels[i], fontsize=16)
#             for j in range(0, len(merger_data)):
#                 a = np.log10(merger_data[j]['galaxy_m'+str(m)])
#                 b = np.log10(merger_data[j][str(names[i])+str(m)])
#                 if merger_labels[j]=='Non merger':
#                     x,y,ysig = myrunningmedian(a,b,20)
#                     #print(titles[m],np.average(ysig)/np.sqrt(len(ysig)))
#                     axes['redbin'+str(m)].plot(x, y, color = colours[j], label=merger_labels[j])
#                     axes['redbin'+str(m)].fill_between(x, y-ysig, y+ysig, facecolor=colours[j], alpha=0.25)
#                     if names[i]=='burst_ssfr' and titles[m]==r'$2 < z < 2.5$':
#                         ax2.plot(x, y, color = colours[j], label=merger_labels[j])
#                         ax2.fill_between(x, y-ysig, y+ysig, facecolor=colours[j], alpha=0.25)
#                 else:
#                     if names[i]=='burst_ssfr':
#                         for k in range(0, len(b)):
#                             if b[k]<-10.5 and titles[m]==r'$1 < z < 1.5$':
#                                 b[k] = -10.5
#                                 axes['redbin'+str(m)].arrow(a[k],b[k],0,-0.3, color=colours[j],
#                                                                 head_width=0.015, width=0.001,
#                                                                 head_length=0.1)
#                             elif b[k]<-12:
#                                 b[k] = -12
#                                 axes['redbin'+str(m)].arrow(a[k],b[k],0,-0.3, color=colours[j],
#                                                                 head_width=0.015, width=0.001,
#                                                                 head_length=0.1)
#                     elif names[i]=='sfe_gal' and titles[m]==r'$1 < z < 1.5$':
#                         for k in range(0, len(b)):
#                             if b[k]<-10.2:
#                                 b[k] = -10.2
#                                 axes['redbin'+str(m)].arrow(a[k],b[k],0,-0.3, color=colours[j],
#                                                                 head_width=0.015, width=0.001,
#                                                                 head_length=0.1)
#                     elif names[i]=='gas_frac' and titles[m]==r'$1 < z < 1.5$':
#                         for k in range(0, len(b)):
#                             if b[k]<-2:
#                                 b[k]=-2
#                                 axes['redbin'+str(m)].arrow(a[k],b[k],0,-0.3, color=colours[j],
#                                                                 head_width=0.015, width=0.001,
#                                                                 head_length=0.1)
#                     x,y,ysig = myrunningmedian(a,b,15)
#                     axes['redbin'+str(m)].scatter(a,b, color=colours[j], label=merger_labels[j], marker='.')
#                     axes['redbin'+str(m)].plot(x, y, color = colours[j], linewidth=2.5)
#                     if names[i]=='burst_ssfr' and titles[m]==r'$2 < z < 2.5$':
#                         ax2.scatter(a,b, color=colours[j], label=merger_labels[j], marker='.')
#                         ax2.plot(x, y, color = colours[j], linewidth=2.5)
#             axes['redbin'+str(m)].text(0.05, 0.05, titles[m], transform=axes['redbin'+str(m)].transAxes, fontsize=14,
#             verticalalignment='bottom', bbox=props)
#             axes['redbin'+str(m)].margins(.2)
#             axes['redbin'+str(m)].set_xlim([9.3,11.9])
#
#         axes['redbin'+str(len(titles)-1)].set_xlabel(r'$\log(M_{*})$', fontsize=16)
#
#         axes['redbin0'].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                ncol=3, mode="expand", borderaxespad=0., prop={'size': 13})
#         fig.tight_layout()
#         fig.savefig('merger_'+str(names[i])+'_h2.png', format='png', dpi=200)
#     ax2.set_xlabel(r'$\log(M_{*})$', fontsize=16)
#     ax2.set_ylabel(r'$\log(sSFR)$', fontsize=16)
#     ax2.set_ylim([-10.0, -8.0])
#     ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
#                 ncol=3, mode="expand", borderaxespad=0., prop={'size': 13})
#     fig2.tight_layout()
#     fig2.savefig(str(results_folder)+'merger_ssfr_forreport.eps', format='eps', dpi=250)

# def merger_fractionPlot(redshift_mer, redshift_all):
#     titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
#     markers = ['o','v', 's']
#     fig = plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
#     ax = fig.add_subplot(1,1,1)
#     ax.set_xlabel(r'$\log(M_{*})$', fontsize=16)
#     ax.set_ylabel(r'Fraction of mergers (%)', fontsize=16)
#     n_mbins = 10
#     for m in range(0, len(titles)):
#         print(len(redshift_mer['galaxy_m'+str(m)]))
#         print(len(redshift_all['galaxy_m'+str(m)]))
#         m_mergers = np.log10(redshift_mer['galaxy_m'+str(m)])
#         m_all = np.log10(redshift_all['galaxy_m'+str(m)])
#         maxs = np.array([m_mergers.max(),m_all.max()])
#         mins = np.array([m_mergers.min(),m_all.min()])
#         mbins = np.linspace(mins.min(), maxs.max(), n_mbins)
#         delta = mbins[1] - mbins[0]
#
#         digi_mer = np.digitize(m_mergers, mbins, right=True)
#         digi_all = np.digitize(m_all, mbins, right=True)
#         binco_mer = np.bincount(digi_mer)
#         binco_all = np.bincount(digi_all)
#         bin_cent = mbins - delta/2
#         frac_mer = []
#         bin_frac =[]
#         for i in range(0, len(binco_mer)):
#             if binco_mer[i] != 0 and binco_all[i] != 0:
#                 frac_mer.append(100*binco_mer[i]/(binco_mer[i]+binco_all[i]))
#                 bin_frac.append(bin_cent[i])
#         ax.plot(bin_frac, frac_mer, label=titles[m], linestyle='--', marker=markers[m])
#     ax.legend(loc='best', prop={'size': 12})
#     fig.tight_layout()
#     fig.savefig(str(results_folder)+'fraction_mergers.png', format='png', dpi=200)

def statsMergers(mergers, sf_galaxies, nbins, printresults = True, plot=False):
    ylabels = [r'$\log(sSFR)$',r'$\log(F_{H_2})$',r'$\log(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    stats_results = {}
    for m in range(0, len(names)):
        stats_results[names[m]] = {}
        for i in range(0, len(titles)):
            stats_results[names[m]][titles[i]] = {}
            aft_m = np.log10(np.asarray([i.m_gal[2] for i in mergers]))
            bef_m = np.log10(np.asarray([i.m_gal[0] for i in mergers]))
            msq = np.log10(np.asarray([i.m_gal for i in sf_galaxies]))
            if m==0:
                aft = np.log10(np.asarray([i.sfr_gal[2] for i in mergers]))
                bef = np.log10(np.asarray([i.sfr_gal[0] for i in mergers]))
                msq_m = np.log10(np.asarray([i.ssfr_gal for i in sf_galaxies]))
            elif m==1:
                aft = np.log10(np.asarray([i.fgas_gal[2] for i in mergers]))
                bef = np.log10(np.asarray([i.fgas_gal[0] for i in mergers]))
                msq_m = np.log10(np.asarray([i.fgas_gal for i in sf_galaxies]))
            elif m==2:
                aft = np.log10(np.asarray([i.sfe_gal[2] for i in mergers]))
                bef = np.log10(np.asarray([i.sfe_gal[0] for i in mergers]))
                msq_m = np.log10(np.asarray([i.sfe_gal for i in sf_galaxies]))
            maxs = np.array([aft_m.max(),bef_m.max(),msq_m.max()])
            mins = np.array([aft_m.min(),bef_m.min(),msq_m.min()])
            bins = np.linspace(mins.min(), maxs.max(), nbins)
            delta = bins[1] - bins[0]
            bin_cent = bins - delta/2
            print(stats_results[names[m]])
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
            for i in range(0, len(titles)):
                print('----------------------------------')
                print('Redshift bin considered: '+str(titles[i]))
                for j in range(0, nbins-1):
                    print('................')
                    print('Mass bin center: '+str(stats_results[names[m]][titles[i]]['bin_cent'][j]))
                    print('p-value from KS 2-sample test: '+str(stats_results[names[m]][titles[i]]['merger_pvalue'][j]))
            print('Statistical significance of difference between after and before merger data')
            for i in range(0, len(titles)):
                print('----------------------------------')
                print('Redshift bin considered: '+str(titles[i]))
                for j in range(0, nbins-1):
                    print('................')
                    print('Mass bin center: '+str(stats_results[names[m]][titles[i]]['bin_cent'][j]))
                    print('p-value from KS 2-sample test: '+str(stats_results[names[m]][titles[i]]['aftvsbef_pvalue'][j]))

def distanceMSQ(merger_data, nbins):
    ylabels = [r'$\Delta_{MSQ}(sSFR)$',r'$\Delta_{MSQ}(f_{H_2})$',r'$\Delta_{MSQ}(SFE)$']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    merger_labels = ['After merger','Before merger','Non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
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
            aft = np.log10(merger_data[0][str(names[i])+str(m)])
            aft_m = np.log10(merger_data[0]['galaxy_m'+str(m)])
            bef = np.log10(merger_data[1][str(names[i])+str(m)])
            bef_m = np.log10(merger_data[1]['galaxy_m'+str(m)])
            msq_m = np.log10(merger_data[2]['galaxy_m'+str(m)])
            msq = np.log10(merger_data[2][str(names[i])+str(m)])
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

statsMergers(mergers, sf_galaxies, 5)
#distanceMSQ(merger_data, 10)
after_before_vs_msqPlots(mergers, sf_galaxies)
#merger_fractionPlot(redshift_mer, redshift_all)
