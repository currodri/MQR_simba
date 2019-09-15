#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

This code is an example analysis of the results from the mergerFinder and quenchingFinder code. In this case, the 
analysis performed is the study of the merger population with respect to the star-forming population.

The version here detailed provides the merger statistics plots given in Rodriguez et al. (2019).
@author: currorodriguez
"""

# Import required libraries
import numpy as np
import random
from decimal import Decimal
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns
sns.set(style="ticks")
import cPickle as pickle
import sys

MODEL = sys.argv[1]  # e.g. m50n512
WIND = sys.argv[2]  # e.g. s50 for Simba

# Import other codes
from galaxy_class import GalaxyData, Merger
results_folder = '../mergers/%s/' % (MODEL) # You can change this to the folder where you want your resulting plots
data_file = '/home/curro/quenchingSIMBA/code/SH_Project/mandq_results_%s.pkl' % (MODEL) # File holding the mergerFinder and quenchingFinder info of galaxies

# Extract data from mergers and quenching pickle files
print('Loading pickle file with data...')
obj = open(data_file, 'rb')
data = pickle.load(obj)
obj.close()
galaxies = data['galaxies']
max_redshift_mergers = data['max_redshift_mergers']
print('Data extracted from pickle file!')

def lsfr_condition(type, galaxy, i, d_indx):
    if d_indx != None:
        if type == 'start':
            lsfr = np.log10(1/(galaxy.t[d_indx][i]))-9
        elif type == 'end':
            lsfr  = np.log10(0.2/(galaxy.t[d_indx][i]))-9
    else:
        lsfr = 0
    return lsfr

def plotmedian(x,y,yflag=[],c='k',ltype='--',lw=3,stat='median',bins=8,label=None,pos=None,boxsize=-1):
    if len(yflag) != len(x):
        #print 'Plotmedian: No flag provided, using all values'
        xp = x
        yp = y
    else:
        xp = x[yflag]
        yp = y[yflag]
    b = bins
    # bins<0 sets bins such that there are equal numbers per bin
    if not isinstance(bins, np.ndarray):
        bins = np.arange(0.999*min(xp),1.001*max(xp),(max(xp)-min(xp))/(bins))
    else:
        b = len(bins)-1
    bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bins,statistic=stat)
    bin_std, useless1, useless2 = stats.binned_statistic(xp,yp,bins=bins,statistic='std')
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])

    if boxsize > 0:  # determine cosmic variance over 8 octants, plot errorbars
        if len(yflag) != len(x): posp = pos
        else: posp = pos[yflag]
        pos = np.floor(posp/(0.5*boxsize)).astype(np.int)
        gal_index = pos[:,0] + pos[:,1]*2 + pos[:,2]*4
        bin_oct = np.zeros((abs(b),8))
        for i0 in xrange(8):
            xq = xp[gal_index==i0]
            yq = yp[gal_index==i0]
            bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_edges,statistic=stat)
        bin_oct  = np.ma.masked_invalid(bin_oct)
        var  = np.ma.std(bin_oct, axis=1)
    elif boxsize == -1:
        var = []
        for i0 in range(len(bin_edges)-1):
            xq = xp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
            yq = yp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
    else:
        var = np.zeros(len(bin_cent))
    return bin_cent,bin_means,var,bin_std

def compare_MergMSQ(galaxies, nbins):
    ylabels = [r'$\log$(sSFR[yr$^{-1}$])',r'$\log(f_{H_2})$',r'$\log$(SFE[yr$^{-1}$])']
    names = ['burst_ssfr','gas_frac','sfe_gal']
    merger_labels = ['Merger','MSQ non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    colours = ['y','k']
    colour_lines = ['r']
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    d = {}
    d['bhm'] = []
    d['pos'] = {}
    for label in ylabels:
        d[label] = {}
        for mlabel in merger_labels:
            d[label][mlabel] = []
            d[mlabel] = []
            for i in range(0, len(titles)):
                d[label][mlabel].append([])
                d[mlabel].append([])       
    for mlabel in merger_labels:
        d['pos'][mlabel] = []
        for i in range(0, len(titles)):
            d['pos'][mlabel].append([])
    for i in range(0, len(titles)):
        d['bhm'].append([])
    for gal in galaxies:
        if gal.mergers:
            mergs = np.asarray([int(merg.indx) for merg in gal.mergers]) + 1
        else:
            mergs = np.array([])
        for k in range(0, len(gal.z)):
            zpos = None
            mpos = None
            for i in range(0, len(titles)):
                if zlimits[i][0] <= gal.z[k] < zlimits[i][1]:
                    zpos = i
            if np.any(mergs==k) and gal.bh_m[k]>0:
                mpos = 0
            elif (gal.sfr[0][k]/gal.m[0][k]) >= (10**lsfr_condition('end',gal,k,0)) and gal.h2_gas[k]>0:
                mpos = 1
            if zpos != None and mpos != None:
                d[ylabels[0]][merger_labels[mpos]][zpos].append(np.log10(gal.sfr[0][k]/gal.m[0][k]))
                d[ylabels[1]][merger_labels[mpos]][zpos].append(np.log10(gal.h2_gas[k]/gal.m[0][k]))
                d[ylabels[2]][merger_labels[mpos]][zpos].append(np.log10(gal.sfr[0][k]/gal.h2_gas[k]))
                bhm = np.log10(gal.bh_m[k]/gal.m[0][k])
                if bhm >= -2.0:
                    bhm = -2.0
                elif bhm <= -4.0:
                    bhm = -4.0
                if mpos==0:
                    d['bhm'][zpos].append(bhm)
                d[merger_labels[mpos]][zpos].append(np.log10(gal.m[0][k]))
                d['pos'][merger_labels[mpos]][zpos].append(gal.pos[k])
    
    fig2, axes2 = plt.subplots(len(ylabels), 1, sharex=True, num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
    fig2.subplots_adjust(hspace=0)
    colours2 = ['b','r','tab:orange']
    ylabels2 = [r'$\Delta_{MSQ}$(sSFR[yr$^{-1}$])',r'$\Delta_{MSQ}(f_{H_2})$',r'$\Delta_{MSQ}$(SFE[yr$^{-1}$])']
    axes2[len(ylabels2)-1].set_xlabel(r'$\log(M_{*}[M_{\odot}])$', fontsize=16)
    for l in range(0,len(ylabels)):
        fig, axes = plt.subplots(len(titles), 1, sharex=True, num=None, figsize=(8,11), dpi=80, facecolor='w', edgecolor='k',)
        axes2[l].set_ylabel(ylabels2[l], fontsize=16)
        axes2[l].plot([9.5,12.0],[0.0,0.0], 'k--')
        axes2[l].tick_params(labelsize=12)
        axes2[l].set_ylim([-0.3,0.85])
        axes2[l].set_xlim([9.5,11.6])
        for i in range(0, len(titles)):
            axes[i].set_ylabel(ylabels[l], fontsize=16)
            axes[i].tick_params(labelsize=12)
            mer = np.asarray(d[ylabels[l]][merger_labels[0]][i])
            mer_m = np.asarray(d[merger_labels[0]][i])
            mer_pos = np.asarray(d['pos'][merger_labels[0]][i])
            msq = np.asarray(d[ylabels[l]][merger_labels[1]][i])
            msq_m =  np.asarray(d[merger_labels[1]][i])
            msq_pos = np.asarray(d['pos'][merger_labels[1]][i])
            lbh = np.asarray(d['bhm'][i])
            sc = axes[i].scatter(mer_m,mer, c=lbh,cmap='plasma', label=merger_labels[0],
                    marker='.', s=30.0, alpha=0.7)

            bins = np.arange(0.999*min(mer_m),11.5,(11.5-min(mer_m))/(nbins))
            bins = np.concatenate((bins,np.array([12.0])))
            mer_cen,mer_median,mer_var,mer_std = plotmedian(mer_m,mer,bins=bins,pos=mer_pos,boxsize=data['boxsize_in_kpccm'])
            msq_cen,msq_median,msq_var,msq_std = plotmedian(msq_m,msq,bins=bins,pos=msq_pos,boxsize=data['boxsize_in_kpccm'])

            axes[i].plot(mer_cen, mer_median, color = colour_lines[0], linewidth=2.5)
            axes[i].plot(msq_cen, msq_median, color=colours[1])
            axes[i].fill_between(msq_cen, msq_median-msq_std, msq_median+msq_std, facecolor=colours[1], alpha=0.25)

            axes[i].text(0.05, 0.05, titles[i], transform=axes[i].transAxes, fontsize=14,
                            verticalalignment='bottom', bbox=props)
            axes[i].margins(.2)
            axes[i].set_xlim([9.5,12.0])

            distance = mer_median - msq_median
            distance_std = np.sqrt(mer_var**2+msq_var**2)
            axes2[l].plot(mer_cen, distance, label=titles[i], color=colours2[i])
            axes2[l].fill_between(mer_cen, distance-distance_std, distance+distance_std, facecolor=colours2[i], alpha=0.25)

        fig.subplots_adjust(hspace=0)
        cb = fig.colorbar(sc, ax=axes.ravel().tolist(), orientation='horizontal', pad=0.08)
        cb.set_label(label=r'$\log(M_{BH}/M_*)$', fontsize=16)
        axes[len(titles)-1].set_xlabel(r'$\log(M_{*}[M_{\odot}])$', fontsize=16)
        fig.savefig(str(results_folder)+'merger_'+str(names[l])+'.png', format='png', dpi=200, bbox_inches='tight')
    axes2[2].set_xlabel(r'$\log(M_{*}[M_{\odot}])$', fontsize=16)
    axes2[1].legend(loc='best', prop={'size': 12})
    fig2.savefig(str(results_folder)+'distance_msq.png', format='png', dpi=200, bbox_inches='tight')


def compare_MergMSQ2(galaxies, nbins):
    ylabels = [r'$\log$(sSFR[yr$^{-1}$])',r'$\log(M_{H_2}/M_{HI})$',r'$\log$(SFE[yr$^{-1}$])']
    names = ['burst_ssfr','gas_ratio','sfe_gal']
    merger_labels = ['Merger','MSQ non merger']
    titles = [r'$0 < z < 0.5$',r'$1 < z < 1.5$',r'$2 < z < 2.5$']
    zlimits = [[0.0, 0.5], [1.0, 1.5], [2.0, 2.5]]
    colours = ['y','k']
    colour_lines = ['r']
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='k')
    d = {}
    d['bhm'] = []
    d['pos'] = {}
    for label in ylabels:
        d[label] = {}
        for mlabel in merger_labels:
            d[label][mlabel] = []
            d[mlabel] = []
            for i in range(0, len(titles)):
                d[label][mlabel].append([])
                d[mlabel].append([])       
    for mlabel in merger_labels:
        d['pos'][mlabel] = []
        for i in range(0, len(titles)):
            d['pos'][mlabel].append([])
    for i in range(0, len(titles)):
        d['bhm'].append([])
    for gal in galaxies:
        if gal.mergers:
            mergs = np.asarray([int(merg.indx) for merg in gal.mergers]) + 1
        else:
            mergs = np.array([])
        for k in range(0, len(gal.z)):
            zpos = None
            mpos = None
            for i in range(0, len(titles)):
                if zlimits[i][0] <= gal.z[k] < zlimits[i][1]:
                    zpos = i
            if np.any(mergs==k):#and gal.bh_m[k]>0:
                mpos = 0
            elif (gal.sfr[0][k]/gal.m[0][k]) >= (10**lsfr_condition('end',gal,k,0)) and gal.h2_gas[k]>0 and gal.h1_gas[k]>0:
                mpos = 1
            if zpos != None and mpos != None:
                d[ylabels[0]][merger_labels[mpos]][zpos].append(np.log10(gal.sfr[0][k]/gal.m[0][k]))
                d[ylabels[1]][merger_labels[mpos]][zpos].append(np.log10(gal.h2_gas[k]/gal.h1_gas[k]))
                d[ylabels[2]][merger_labels[mpos]][zpos].append(np.log10(gal.sfr[0][k]/gal.h2_gas[k]))
                bhm = np.log10(gal.bh_m[k]/gal.m[0][k])
                if bhm >= -2.0:
                    bhm = -2.0
                elif bhm <= -4.0:
                    bhm = -4.0
                if mpos==0:
                    d['bhm'][zpos].append(bhm)
                d[merger_labels[mpos]][zpos].append(np.log10(gal.m[0][k]))
                d['pos'][merger_labels[mpos]][zpos].append(gal.pos[k])
    
    fig2, axes2 = plt.subplots(len(ylabels), 1, sharex=True, num=None, figsize=(8, 9), dpi=80, facecolor='w', edgecolor='k')
    fig2.subplots_adjust(hspace=0)
    colours2 = ['b','r','tab:orange']
    ylabels2 = [r'$\Delta_{MSQ}$(sSFR[yr$^{-1}$])',r'$\Delta_{MSQ}(M_{H_2}/M_{HI})$',r'$\Delta_{MSQ}$(SFE[yr$^{-1}$])']
    axes2[len(ylabels2)-1].set_xlabel(r'$\log(M_{*}[M_{\odot}])$', fontsize=16)
    for l in range(0,len(ylabels)):
        fig, axes = plt.subplots(len(titles), 1, sharex=True, num=None, figsize=(8,11), dpi=80, facecolor='w', edgecolor='k',)
        axes2[l].set_ylabel(ylabels2[l], fontsize=16)
        axes2[l].plot([9.5,12.0],[0.0,0.0], 'k--')
        axes2[l].tick_params(labelsize=12)
        axes2[l].set_ylim([-0.3,0.85])
        axes2[l].set_xlim([9.5,11.6])
        for i in range(0, len(titles)):
            axes[i].set_ylabel(ylabels[l], fontsize=16)
            axes[i].tick_params(labelsize=12)
            mer = np.asarray(d[ylabels[l]][merger_labels[0]][i])
            mer_m = np.asarray(d[merger_labels[0]][i])
            mer_pos = np.asarray(d['pos'][merger_labels[0]][i])
            msq = np.asarray(d[ylabels[l]][merger_labels[1]][i])
            msq_m =  np.asarray(d[merger_labels[1]][i])
            msq_pos = np.asarray(d['pos'][merger_labels[1]][i])
            lbh = np.asarray(d['bhm'][i])
            sc = axes[i].scatter(mer_m,mer, c=lbh,cmap='plasma', label=merger_labels[0],
                    marker='.', s=30.0, alpha=0.7)

            bins = np.arange(0.999*min(mer_m),11.5,(11.5-min(mer_m))/(nbins))
            bins = np.concatenate((bins,np.array([12.0])))
            mer_cen,mer_median,mer_var,mer_std = plotmedian(mer_m,mer,bins=bins,pos=mer_pos,boxsize=data['boxsize_in_kpccm'])
            msq_cen,msq_median,msq_var,msq_std = plotmedian(msq_m,msq,bins=bins,pos=msq_pos,boxsize=data['boxsize_in_kpccm'])

            axes[i].plot(mer_cen, mer_median, color = colour_lines[0], linewidth=2.5)
            axes[i].plot(msq_cen, msq_median, color=colours[1])
            axes[i].fill_between(msq_cen, msq_median-msq_std, msq_median+msq_std, facecolor=colours[1], alpha=0.25)

            axes[i].text(0.05, 0.05, titles[i], transform=axes[i].transAxes, fontsize=14,
                            verticalalignment='bottom', bbox=props)
            axes[i].margins(.2)
            axes[i].set_xlim([9.5,12.0])

            distance = mer_median - msq_median
            distance_std = np.sqrt(mer_var**2+msq_var**2)
            axes2[l].plot(mer_cen, distance, label=titles[i], color=colours2[i])
            axes2[l].fill_between(mer_cen, distance-distance_std, distance+distance_std, facecolor=colours2[i], alpha=0.25)

        fig.subplots_adjust(hspace=0)
        cb = fig.colorbar(sc, ax=axes.ravel().tolist(), orientation='horizontal', pad=0.08)
        cb.set_label(label=r'$\log(M_{BH}/M_*)$', fontsize=16)
        axes[len(titles)-1].set_xlabel(r'$\log(M_{*}[M_{\odot}])$', fontsize=16)
        fig.savefig(str(results_folder)+'merger_'+str(names[l])+'.png', format='png', dpi=200, bbox_inches='tight')
    axes2[1].legend(loc='best', prop={'size': 12})
    fig2.savefig(str(results_folder)+'distance_msq.png', format='png', dpi=200, bbox_inches='tight')

compare_MergMSQ(galaxies,10)