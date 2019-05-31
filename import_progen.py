#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 December 2018

@author: currorodriguez
"""

# Import required libraries
import numpy as np
import sys
def numGal(filename):
    data0 = []
    delimiters = []
    i = 0

    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    ngal = len(delimiters)

    return(ngal)

def importApp(foldername):
    sys.path.insert(0, foldername)
    # Generate initial numpy array from txt file
    filename = '/progen_mass.txt'
    data0 = []
    delimiters = []
    i = 0

    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    d = {} # Create dictionary to store the diferent galaxies variables

    ngal = len(delimiters)

    for j in range(0,ngal):
        if j==0:
            d["galaxy_t" + str(j)] = np.zeros(delimiters[j])
            d["m_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["galaxy_t" + str(j)][k] = float(data0[k][0])
                d["m_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["galaxy_t" + str(j)] = np.zeros(n-1)
            d["m_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["galaxy_t" + str(j)][k] = float(data0[delimiters[j-1]+1+k][0])
                d["m_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])
    # Generate initial numpy array from txt file
    filename = '/progen_sfr.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1


    for j in range(0,ngal):
        if j==0:
            d["sfr_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["sfr_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["sfr_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["sfr_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])

    # Generate initial numpy array from txt file
    filename = '/progen_fgas.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(0,ngal):
        if j==0:
            d["fgas_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["fgas_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["fgas_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["fgas_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])

    # Generate initial numpy array from txt file
    filename = '/progen_h1_gal.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(0,ngal):
        if j==0:
            d["h1_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["h1_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["h1_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["h1_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])
    # Generate initial numpy array from txt file
    filename = '/progen_h2_gal.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(0,ngal):
        if j==0:
            d["h2_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["h2_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["h2_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["h2_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])

    # Generate initial numpy array from txt file
    filename = '/progen_z.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(0,ngal):
        if j==0:
            d["z_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["z_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["z_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["z_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])

    # Generate initial numpy array from txt file
    filename = '/progen_sfe.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(0,ngal):
        if j==0:
            d["sfe_gal" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["sfe_gal" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["sfe_gal" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["sfe_gal" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])

# Generate initial numpy array from txt file
    filename = '/progen_gal_type.txt'
    data0 = []
    delimiters = []
    i = 0
    j = 0
    row = 0
    n = 0
    k = 0
    with open(filename) as readfile:
        for row in readfile:
            a = row.split()
            data0.append(a)
            if a==['#']: # Each galaxy data is separated using #
                delimiters.append(i)
            i = i+1

    for j in range(ngal):
        if j==0:
            d["gal_type" + str(j)] = np.zeros(delimiters[j])
            for k in range(0,delimiters[j]):
                d["gal_type" + str(j)][k] = float(data0[k][1])
        else:
            n = delimiters[j]-delimiters[j-1]
            d["gal_type" + str(j)] = np.zeros(n-1)
            for k in range(0, n-1):
                d["gal_type" + str(j)][k] = float(data0[delimiters[j-1]+1+k][1])
    return (d, ngal)
print('Data extracted from txt files.')
