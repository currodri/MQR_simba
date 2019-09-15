#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 1 10:22:59 2019

@author: Curro Rodriguez Montero, School of Physics and Astronomy,
            University of Edinburgh, JCMB, King's Buildings

For questions about the code:
s1650043@ed.ac.uk
"""
"""Import some necessary packages"""
import numpy as np

""" Define classes """

class GalaxyData:
    def __init__(self,progen_id,sfr,m,z,t,h1_gas,h2_gas,bh_m,bhar,local_den,g_type,pos,caesar_id):
        self.interpolation = False
        self.progen_id = int(progen_id)
        self.sfr = [sfr,0]
        self.ssfr = [0,0]
        self.m = [m,0]
        self.z = z
        self.t = [t,0]
        self.h1_gas = h1_gas
        self.h2_gas = h2_gas
        self.bh_m = bh_m
        self.bhar = bhar
        self.local_den = local_den
        self.g_type = g_type
        self.pos = pos
        self.caesar_id = caesar_id
        self.mergers = []
        self.quenching = []
        self.rejuvenations = []
        self.mags = []
        self.scs = []
    def get_ssfr(self):
        if self.interpolation:
            self.ssfr[1] = self.sfr[1]/self.m[1]
        else:
            self.ssfr[0] = self.sfr[0]/self.m[0]
    def get_fgas(self):
        self.fgas = self.h2_gas/self.m
    def get_sfe(self):
        self.sfe = self.sfr/self.h2_gas
    def interpolated_data(self,sfr_new,m_new,t_new):
        self.sfr[1] = np.asarray(sfr_new)
        self.m[1] = np.asarray(m_new)
        self.t[1] = np.asarray(t_new)

class Quench:
    def __init__(self, above9):
        self.above9 = above9
        self.below11 = None
        self.quench_time = None
        self.indx = None

class Merger:
    def __init__(self,indx,merger_ratio,fgas_boost):
        self.indx = indx
        self.merger_ratio = merger_ratio
        self.fgas_boost = fgas_boost

class Magnitude:
    def __init__(self):
        self.filtername = None
        self.wave_eff = None
        self.z = []
        self.Abs = []
        self.App = []

class SuperColour:
    def __init__(self):
        self.sc_number = None
        self.z = []
        self.values = []