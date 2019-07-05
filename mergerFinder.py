#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: Curro Rodriguez Montero, School of Physics and Astronomy,
            University of Edinburgh, JCMB, King's Buildings

For questions about the code:
s1650043@ed.ac.uk
"""

"""Import some necessary packages"""
import numpy as np
import pylab as plt
from scipy import stats
from quenchingFinder import sfr_condition_2, GalaxyData

"""Classes defined"""
class Merger:
    def __init__(self,id, sfr_gal, sfe_gal, z_gal, galaxy_t, m_gal, fgas_gal, gal_type, gal_pos, merger_ratio, fgas_boost):
        self.sfr_gal = sfr_gal
        self.ssfr_gal = (sfr_gal/m_gal)
        self.sfe_gal = sfe_gal
        self.z_gal = z_gal
        self.galaxy_t = galaxy_t
        self.m_gal = m_gal
        self.fgas_gal = fgas_gal
        self.id = id
        self.type = gal_type
        self.gal_pos = gal_pos
        self.merger_ratio = merger_ratio
        self.fgas_boost = fgas_boost
###########################################################################################
"""
FUNCTION THAT DEFINES THE CONDITIONS FOLLOWED TO DETECT A MERGER
"""
def merger_condition(mass_list, index, merger_ratio, mass_limit):
    condition = False
    diff = (mass_list[index+1]-mass_list[index])/mass_list[index]
    diff2 = abs((mass_list[index+2]-mass_list[index])/mass_list[index])
    diff3 = abs((mass_list[index+1]-mass_list[index-1])/mass_list[index-1])
    diff4 = abs((mass_list[index+3]-mass_list[index])/mass_list[index])
    if diff>=merger_ratio and diff2>=merger_ratio and diff-diff3 < 0.001 and diff4>=merger_ratio and mass_list[index]>=mass_limit:
        condition = True
    return (condition, diff)
###########################################################################################
"""
MAIN FUNCTION FOR THE MERGER FINDER OVER THE GALAXIES

ARGUMENTS

galaxies ======= dictionary containing all the galaxies with their properties
                    (ssfr in this case should be sfr instead)
merger_ratio === the ratio above which the code looks for mergers, i.e. R=4:1 would be merger_ratio=0.2
mass_limit ===== minimum mass of final galaxy at which the code looks for mergers
redshift_limit = maximum redshift at which the code looks for mergers

"""
def merger_finder(galaxies, merger_ratio, mass_limit, redshift_limit):
    mergers = []
    sf_galaxies = []
    for galaxy in range(0, len(galaxies)):
        gal = galaxies[galaxy]
        mass = gal.m_gal
        sfr = gal.sfr_gal
        fgas = gal.fgas_gal
        z = gal.z_gal
        sfe = gal.sfe_gal
        id = gal.id
        type = gal.type
        time = gal.galaxy_t
        pos = gal.gal_pos
        for i in range(1, len(mass)-3):
            if z[i]<=redshift_limit:
                condition,ratio = merger_condition(mass, i, merger_ratio, mass_limit)
                sfcondition = sfr_condition_2('end', gal, i)
                ssfr = sfr[i]/mass[i]
                if condition == True and ssfr>=(10**sfcondition) and fgas[i]>0:
                    boost = (fgas[i+1]-fgas[i-1])/fgas[i-1]
                    # Save data at the merger, after and before
                    merger = Merger(id,sfr[i-1:i+2],sfe[i-1:i+2],z[i-1:i+2],time[i-1:i+2],mass[i-1:i+2],fgas[i-1:i+2],type[i-1:i+2],pos[i-1:i+2],ratio, boost)
                    mergers.append(merger)
                else:
                    if mass[i]>mass_limit and ssfr>=(10**sfcondition):
                        sf_gal = GalaxyData(id,sfr[i],sfe[i],z[i],time[i],mass[i],fgas[i],type[i], pos[i])
                        # Add star forming galaxy to the list
                        sf_galaxies.append(sf_gal)
    print('Star-forming main sequence and mergers found up to z = '+str(redshift_limit))
    return(mergers, sf_galaxies)


##########################################################################################
"""
EXTRA FUNCTIONS USEFUL FOR THE ANALYSIS OF THE RESULTS
"""

def histedges_equalN(x, nbin):
    npt = len(x)
    return np.interp(np.linspace(0, npt, nbin + 1),
                     np.arange(npt),
                     np.sort(x))

def plotmedian(x,y,yflag=[],c='k',ltype='--',lw=3,stat='median',ax='plt',bins=8,label=None,pos=None,boxsize=-1, bin_choosen=0):
    if len(yflag) != len(x):
        #print 'Plotmedian: No flag provided, using all values'
        xp = x
        yp = y
    else:
        xp = x[yflag]
        yp = y[yflag]
    # bins<0 sets bins such that there are equal numbers per bin
    #if bins < 0: bin_edges = histedges_equalN(xp,-bins)
    #else: bin_edges = np.arange(0.999*min(xp),1.001*max(xp),(max(xp)-min(xp))/(bins))
    #else: bin_edges = bin_choosen
    #if bins < 0:	bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bin_edges,statistic=stat)
    #else: bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bins,statistic=stat)
    #bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    #ax.plot(bin_cent, bin_means, ltype, lw=lw, color=c, label=label)
    bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bin_choosen,statistic=stat)
    print(bin_means, bin_edges, binnumber)
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])

    if boxsize > 0:  # determine cosmic variance over 8 octants, plot errorbars
	if len(yflag) != len(x): posp = pos
	else: posp = pos[yflag]
        pos = np.floor(posp/(0.5*boxsize)).astype(np.int)
        gal_index = pos[:,0] + pos[:,1]*2 + pos[:,2]*4
        bin_oct = np.zeros((abs(bins),8))
        for i0 in xrange(8):
            xq = xp[gal_index==i0]
            yq = yp[gal_index==i0]
	    if bins < 0: bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_edges,statistic=stat)
        #else: bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_choosen,statistic=stat)
        else: bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_edges,statistic=stat)
        bin_oct  = np.ma.masked_invalid(bin_oct)
        var  = np.ma.std(bin_oct, axis=1)
        print(bin_oct[:,i0], bin_edges, binnumber)
        print(var)
        #ax.errorbar(bin_cent, bin_means, yerr=[var,var], fmt='o', linewidth=lw, color=c)
    elif boxsize == -1:
        var = []
        for i0 in range(len(bin_edges)-1):
            xq = xp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
            yq = yp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
            var.append(np.ma.std(yq-np.mean(yq)))
        #print 'bins',bin_cent
        #print 'variance',var
        #ax.errorbar(bin_cent, bin_means, yerr=[var,var], fmt='o', linewidth=lw, color=c)
    else:
        var = np.zeros(len(bin_cent))
    return bin_means,var


def plotmedian2(x,y,yflag=[],c='k',ltype='--',lw=3,stat='median',ax='plt',bins=8,label=None,pos=None,boxsize=-1):
    if len(yflag) != len(x):
        #print 'Plotmedian: No flag provided, using all values'
        xp = x
        yp = y
    else:
        xp = x[yflag]
        yp = y[yflag]
    # bins<0 sets bins such that there are equal numbers per bin
    if bins < 0: bin_edges = histedges_equalN(xp,-bins)
    else: bin_edges = np.arange(0.999*min(xp),1.001*max(xp),(max(xp)-min(xp))/(bins))
    if bins < 0:	bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bin_edges,statistic=stat)
    else: bin_means, bin_edges, binnumber = stats.binned_statistic(xp,yp,bins=bins,statistic=stat)
    print(bin_edges)
    print(np.arange(0.999*min(xp),1.001*max(xp),(max(xp)-min(xp))/(bins)))
    bin_cent = 0.5*(bin_edges[1:]+bin_edges[:-1])
    #ax.plot(bin_cent, bin_means, ltype, lw=lw, color=c, label=label)

    if boxsize > 0:  # determine cosmic variance over 8 octants, plot errorbars
	if len(yflag) != len(x): posp = pos
	else: posp = pos[yflag]
        pos = np.floor(posp/(0.5*boxsize)).astype(np.int)
        gal_index = pos[:,0] + pos[:,1]*2 + pos[:,2]*4
        bin_oct = np.zeros((abs(bins),8))
        for i0 in xrange(8):
            xq = xp[gal_index==i0]
            yq = yp[gal_index==i0]
	    if bins < 0: bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_edges,statistic=stat)
	    else: bin_oct[:,i0], bin_edges, binnumber = stats.binned_statistic(xq,yq,bins=bin_edges,statistic=stat)
        bin_oct  = np.ma.masked_invalid(bin_oct)
        if stat=='median':
            var  = np.ma.std(bin_oct, axis=1)
        elif stat=='mean':
            lens = np.zeros(len(bin_oct))
            for i in range(0, len(bin_oct)):
                lens[i] = len(bin_oct[i])
            var  = np.ma.std(bin_oct, axis=1)/np.ma.sqrt(lens)
        #ax.errorbar(bin_cent, bin_means, yerr=[var,var], fmt='o', linewidth=lw, color=c)
    elif boxsize == -1:
        var = []
        for i0 in range(len(bin_edges)-1):
            xq = xp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
            yq = yp[(xp>bin_edges[i0])&(xp<bin_edges[i0+1])]
            var.append(np.ma.std(yq-np.mean(yq)))
        #print 'bins',bin_cent
        #print 'variance',var
        #ax.errorbar(bin_cent, bin_means, yerr=[var,var], fmt='o', linewidth=lw, color=c)
    else:
        var = np.zeros(len(bin_cent))
    return bin_cent,bin_means,var
# Running median through scatter data
def myrunningmedian(x,y,nbins, sigma=True):
    bins = np.linspace(x.min(), x.max(), nbins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_median = [np.median(y[idx==k]) for k in range(0,nbins)]
    running_median = np.asarray(running_median)
    delitems = []
    ynan = np.isnan(running_median)
    for i in range(0, len(ynan)):
        if ynan[i]:
            delitems.append(i)
    running_median = np.delete(running_median, delitems)
    bins = np.delete(bins, delitems)
    bin_cent = bins - delta/2
    if sigma==True:
        running_std = [y[idx==k].std() for k in range(0,nbins)]
        running_std = np.asarray(running_std)
        running_std = np.delete(running_std, delitems)
        return bin_cent, running_median, running_std
    else:
        return bin_cent, running_median
