#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 15:22:59 2018

@author: Curro Rodriguez Montero, School of Physics and Astronomy,
            University of Edinburgh, JCMB, King's Buildings

For questions about the code:
s1650043@ed.ac.uk
"""
"""Import some necessary packages"""
import numpy as np
from scipy import interpolate


"""Classes defined"""
class GalaxyData:
    def __init__(self,id, ssfr_gal, sfe_gal, z_gal, galaxy_t, m_gal, fgas_gal, gal_type):
        self.ssfr_gal = ssfr_gal
        self.sfe_gal = sfe_gal
        self.z_gal = z_gal
        self.galaxy_t = galaxy_t
        self.m_gal = m_gal
        self.fgas_gal = fgas_gal
        self.quenching = []
        self.id = id
        self.type = gal_type
        self.rate = []

class Quench:
    def __init__(self, above9, type):
        self.above9 = above9
        self.type = type
        self.below11 = None
        self.quench_time = None


###########################################################################################
"""
MAIN FUNCTION FOR THE QUENCHING FINDER OVER THE GALAXIES

ARGUMENTS

galaxies ======= dictionary containing all the galaxies with their properties
sfr_condition == method that will be used for the thresholds in star formation and quenching
mass_limit ===== minimum mass of final galaxy at which the code looks for quenching
interpolation == if set to True, the interpolated data is used for the quenching analysis. If
                    set to nothing it is set to False

"""


def quenchingFinder2(galaxies,sfr_condition, mass_limit, interpolation=False):

    sfr_conditions = [sfr_condition_1, sfr_condition_2]
    sfr_condition = sfr_conditions[int(sfr_condition)]
    interpolation_list_of_list = []
    total_quenched = 0

    for i in range(0, len(galaxies)):

        #Galaxy we are considering
        galaxy = galaxies[i]
        lookup_condition = sfr_condition('end', galaxy, -1)
        if not interpolation:
            m = np.log10(galaxy.m_gal[-1])
        else:
            m = np.log10(galaxy.m_gal)
        if galaxy.ssfr_gal[-1]<(10**lookup_condition) and m>=mass_limit and not interpolation:
            total_quenched = total_quenched + 1
            #State of the search
            state = (0, galaxy.galaxy_t[0], None)
            #The state has 3 elements.
            #The first one indicates the stage we are in (initial, pre_quench or quench)
            #The second one has the inital time of the period we are considering
            #The third one has the time of the moment we found a pre_quench


            #Set the number of snapshots to be observed
            last_snapshot = len(galaxy.galaxy_t)

            #Go over each snapshot and save the new data of the galaxy
            for j in range(0, last_snapshot):
                state = analyseState[state[0]](galaxy,j, state, sfr_condition)
            #Check if the last quenching is a valid one:
            if galaxy.quenching and galaxy.quenching[-1].below11 == None:
                del galaxy.quenching[-1]

            # galaxy_interpolated = ssfr_interpolation(galaxy)
            if galaxy.quenching:
                galaxies_interpolated = ssfr_interpolation(galaxy)
                interpolation_list_of_list.extend(galaxies_interpolated)
        elif interpolation:
            total_quenched = total_quenched + 1
            #State of the search
            state = (0, galaxy.galaxy_t[0], None)
            #The state has 3 elements.
            #The first one indicates the stage we are in (initial, pre_quench or quench)
            #The second one has the inital time of the period we are considering
            #The third one has the time of the moment we found a pre_quench


            #Set the number of snapshots to be observed
            last_snapshot = len(galaxy.galaxy_t)

            #Go over each snapshot and save the new data of the galaxy
            for j in range(0, last_snapshot):
                state = analyseState[state[0]](galaxy,j, state, sfr_condition, interpolation=True)

            #Check if the last quenching is a valid one:
            if galaxy.quenching and galaxy.quenching[-1].below11 == None:
                del galaxy.quenching[-1]
    print ('Total number of quenched galaxies at z=0 : '+str(total_quenched))
    return interpolation_list_of_list


###########################################################################################
"""
FUNCTIONS THAT DEFINE THE DIFFERENT STAGES FOR QUENCHING AND REJUVENATION
"""

def initial(galaxy,j,curr_state, sfr_condition, interpolation=False):
    """We check if the ssfr is higher than threshold... if that's the case, then we are
    ready to look for a quench. """

    ssfr_gal, galaxy_t = galaxy.ssfr_gal[j], galaxy.galaxy_t[j]

    current_lssfr = sfr_condition('start', galaxy, j)

    if ssfr_gal > 10**current_lssfr:
        new_state = (1, galaxy_t, None)
    else:
        new_state = (0, None, None)

    return new_state

def readyToLook (galaxy,j,curr_state, sfr_condition, interpolation=False):
    """We are ready to check if ssfr is below threshold"""
    ssfr_gal = galaxy.ssfr_gal[j]

    current_lssfr = sfr_condition('start', galaxy, j)

    if ssfr_gal <= 10**current_lssfr:
        if not interpolation:
            quench = Quench(j-1, galaxy.type[j-1])
        else:
            quench = Quench(j-1, galaxy.type)
        galaxy.quenching.append(quench)

        new_state = (2, curr_state[1], galaxy.galaxy_t[j])
    else:
        new_state = curr_state

    return new_state


def pre_quench (galaxy,j,curr_state, sfr_condition, interpolation=False):
    """There has been a lssfr <= threshold, now let's check for a quench """
    ssfr_gal = galaxy.ssfr_gal[j]

    current_lssfr = sfr_condition('end', galaxy, j)

    if ssfr_gal < 10**current_lssfr:
        #Retrieve the current quench
        quench = galaxy.quenching[-1]
        #Add the point below11 and the length time of the quench
        quench.below11 = j
        quench.quench_time =abs(curr_state[2] - galaxy.galaxy_t[j])

        #Now we look for rejuvenations
        new_state = (3, curr_state[1], None)
    else:
        new_state = curr_state

    return new_state


def quench (galaxy,j,curr_state, sfr_condition, interpolation=False):
    """We have detected a quench and now we are looking for rejuvenations """
    ssfr_gal, galaxy_t = galaxy.ssfr_gal[j], galaxy.galaxy_t[j]

    current_lssfr = sfr_condition('start', galaxy, j)
    time_min = max(curr_state[1], 0.5)

    if galaxy_t > 1.2*time_min:
        #It has passed enough time since the quench was detected, we can have rejuvenations

        if ssfr_gal > 10**current_lssfr:
            #We have found a sign change
            if not interpolation:
                galaxy.rate.append(galaxy.z_gal[j])
                galaxy.rate.append(galaxy_t)
                galaxy.rate.append(galaxy.m_gal[j])
            new_state = (1, galaxy_t, None)
        else:
            new_state = curr_state
    else:
        #Now if there is a rejuvenation, the quench should be discarded.
        if ssfr_gal > 10**current_lssfr:
            #We have found a sign change. Rollback the changes done to the lists.
            del galaxy.quenching[-1]

            #Go back to state readyToLook.
            new_state = (1, galaxy_t, None)
        else:
            #Explore the next snapshot
            new_state = curr_state

    return new_state

##########################################################################################
"""
FUNCTIONS THAT DEFINE THE DIFFERENT INTERPOLATION METHODS THAT CAN BE USE
"""


def ssfr_interpolation(galaxy):

    new_galaxies = []
    aboves = []
    belows = []
    for quench in galaxy.quenching:
        #For each quenching, interpolate the new values creating a new galaxy
        above, below = quench.above9, (quench.below11 + 1)
        limit = 0
        if above - limit < 0 or below + limit >= len(galaxy.galaxy_t):
            limit = min(len(galaxy.galaxy_t) - below, above)
        aboves.append(above-limit)
        belows.append(below+limit)
    aboves = np.asarray(aboves)
    belows = np.asarray(belows)
    above = aboves.min()
    below = belows.max()
    if len(range(above, below+limit,1))>3:
        #If there are at least three points in the quench, then:
        ssfr_gal_non = [galaxy.ssfr_gal[j] for j in range(above-limit, below+limit,1)]
        galaxy_t_non = [galaxy.galaxy_t[j] for j in range(above-limit, below+limit,1)]

        time_new = np.arange(np.amin(galaxy_t_non), np.amax(galaxy_t_non), 0.001)

        tck = interpolate.splrep(galaxy_t_non,ssfr_gal_non, k=3)
        ssfr_new = interpolate.splev(time_new, tck, der=0)

        new_galaxies.append(GalaxyData(galaxy.id, ssfr_new.tolist(), galaxy.sfe_gal[quench.below11],
                                galaxy.z_gal[quench.below11],time_new.tolist(), galaxy.m_gal[quench.above9],
                                galaxy.fgas_gal[quench.above9], quench.type))
    return new_galaxies


##########################################################################################
"""
FUNCTIONS THAT DEFINE THE DIFFERENT THRESHOLDS FOR STAR FORMING AND QUENCHED GALAXIES
"""

def sfr_condition_1(type, galaxy, j):
    if galaxy.z_gal[j]<=2.0:
        a = 0.3
    else:
        a = 0.0
    if type == 'start':
        lsfr = -9.5 + a*galaxy.z_gal[j]
    elif type == 'end':
        lsfr = -11 + a*galaxy.z_gal[j]
    return lsfr

def sfr_condition_2(type, galaxy, j):
    if type == 'start':
        lsfr = np.log10(1/(galaxy.galaxy_t[j]))-9
    elif type == 'end':
        lsfr  = np.log10(0.2/(galaxy.galaxy_t[j]))-9
    return lsfr

analyseState = {0:initial, 1:readyToLook, 2:pre_quench, 3:quench}

##########################################################################################
"""
EXTRA FUNCTIONS USEFUL FOR THE ANALYSIS OF THE RESULTS
"""
def myrunningmedian(x,y,nbins, sigma=True):
    bins = np.linspace(x.min()*0.9, x.max()*1.1, nbins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_median = [np.median(y[idx==k]) for k in range(0,nbins)]
    running_median = np.asarray(running_median)
    delitems = []
    for i in range(0, len(running_median)):
        if np.isnan(running_median[i]):
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
def rejuvenation_rate_calculator(d, rejuvenation_z, count_galaxy_file):

    # Get number of galaxies per snapshot
    num_gal_snap = np.genfromtxt(count_galaxy_file)
    z = d['z_gal0']
    t = d['galaxy_t0']
    zlim = np.amax(rejuvenation_z)
    zlimind = 0
    zbins = []
    tbins = []
    for i in range(0, len(z), 3):
        if z[i]<zlim:
            zbins.append(z[i])
            tbins.append(t[i])
        elif z[i]==zlim:
            zlimind = i
            break
    z = z[0:zlimind-2]
    zbins[0] = zbins[0] * 0.9
    zbins[-1] = zbins[-1] * 1.1
    bin_cent = np.zeros(len(zbins)-1)
    digi = np.digitize(z, bins=zbins, right=True)
    binco = np.bincount(digi)
    histo = np.zeros(len(zbins)-1)
    deltat = np.zeros(len(zbins)-1)
    binco = np.delete(binco, 0)
    for i in range(0, len(digi)):
        index = digi[i]-1
        histo[index] = histo[index] + num_gal_snap[i]
    histo = histo/binco #Average galaxies per redshift bin
    for j in range(0, len(tbins)-1):
        deltat[j] = tbins[j] - tbins[j+1]
        bin_cent[j] = (zbins[j+1] + zbins[j])/2

    digi2 = np.digitize(rejuvenation_z, bins=zbins, right=True)
    binco2 = np.bincount(digi2)
    binco2 = np.delete(binco2, 0)
    rates = binco2/(deltat*histo)
    bin_cent, rates, rates_sig = myrunningmedian(bin_cent, rates, 20)
    return rates, bin_cent, rates_sig

def quenching_histogram(z_init,galaxies,ngal,min_mass, max_mass,quenching_times, redshifts, n_bins):
    z_bins = np.linspace(0.0, np.amax(redshifts)*1.1, n_bins)
    counts = np.zeros(n_bins-1)
    counts_error = np.zeros(n_bins-1)
    times = np.zeros(n_bins-1)
    times_error = np.zeros(n_bins-1)
    delta = z_bins[1] - z_bins[0]
    z_cent = z_bins - delta/2
    z_cent = np.delete(z_cent, 0)
    counts_init = np.zeros(len(z_init)-1)
    z_init_cent = (z_init[:-1]+z_init[1:])/2
    for i in range(0, len(z_init)-1):
        count_m = 0
        count_nm = 0
        for j in range(0, len(redshifts)):
            if z_init[i]<=redshifts[j]<z_init[i+1]:
                count_m = count_m + 1
        for k in range(0, ngal):
            gal = galaxies[k]
            red = gal.z_gal
            mass = gal.m_gal
            for m in range(0, len(red)):
                if z_init[i]<=red[m]<z_init[i+1] and (10**min_mass)<=mass[m]<(10**max_mass):
                    count_nm = count_nm + 1
        if count_m != 0 and count_nm != 0:
            counts_init[i] = count_m/count_nm

    for i in range(0, n_bins-1):
        n_gals = 0
        t = []
        counts_ave = []
        for j in range(0, len(redshifts)):
            if z_bins[i]<=redshifts[j]<z_bins[i+1]:
                t.append(quenching_times[j])
        for k in range(0, len(counts_init)):
            if z_bins[i]<=z_init_cent[k]<z_bins[i+1]:
                counts_ave.append(counts_init[k])
        counts_ave = np.asarray(counts_ave)
        counts[i] = np.average(counts_ave)
        counts_error[i] = np.std(counts_ave)/(np.sqrt(len(counts_ave)))
        t = np.asarray(t)
        times[i] = np.average(t)
        times_error[i] = np.std(t)/(np.sqrt(len(t)))
    return(z_cent, counts, counts_error, times, times_error)


    # z = d['z_gal0']
    # zlim = np.amax(redshifts)
    # zlimind = 0
    # zbins = []
    # for i in range(0, len(z), 3):
    #     if z[i]<zlim:
    #         zbins.append(z[i])
    #     elif z[i]==zlim:
    #         zlimind = i
    #         break
    # z = z[0:zlimind-2]
    # zbins[0] = zbins[0] * 0.9
    # zbins[-1] = zbins[-1] * 1.1
    # bin_cent = np.zeros(len(zbins)-1)
    # digi = np.digitize(redshifts, bins=zbins, right=True)
    # binco = np.bincount(digi)
    # histo = np.zeros(len(zbins)-1)
    # deltat = np.zeros(len(zbins)-1)
    # binco = np.delete(binco, 0)
    # for i in range(0, len(digi)):
    #     index = digi[i]-1
    #     histo[index] = histo[index] + quenching_times[i]/binco[index]
    # for j in range(0, len(zbins)-1):
    #     bin_cent[j] = (zbins[j+1] + zbins[j])/2
    # return histo, binco, bin_cent