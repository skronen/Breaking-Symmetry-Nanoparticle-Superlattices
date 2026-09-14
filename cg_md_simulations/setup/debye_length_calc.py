#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 23 13:52:04 2025

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import e, epsilon_0, Avogadro, Boltzmann

epsilon = 78.4
T = 298
K = (epsilon_0 * epsilon * Boltzmann * T / 2/ Avogadro / 1000 / e**2)**0.5

debye = lambda molar : K/np.sqrt(molar) * 10**9 #convert to nm

def convert_conc_to_debye(concs):
    concs = np.array(concs)/1000 #convert M to mM
    deb = debye(concs)
    return deb * 10 #convert to Angstroms
    
