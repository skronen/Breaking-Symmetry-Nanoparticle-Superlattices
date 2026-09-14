#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 13:30:41 2024

@author: skronen
"""

import numpy as np

sigma = 1# 1.8 #nm
kuhn_mw = 137 #g/mol
kuhn_mw = sigma/1.1 * kuhn_mw
T = 298 #K



def convert_np_diam(diam):
    np_diam = diam/sigma
    return np_diam

def convert_chain_len(MW):
    # MW is polymer MW in g/mol or Dalton
    l = int(np.round(MW/kuhn_mw, 0))
    return l
    

def convert_gd(diam, gd):
    SA = np.pi * diam**2
    n_chains = int(np.round(gd * SA, 0))
    gdens = gd * sigma**2
    
    return n_chains, gdens

def convert_eng(eng):
    #eng in kJ/mol
    epsilon = T * 8.314462618e-3 #kb in kJ/mol.K
    eps_LJ = eng/epsilon
    return eps_LJ
    

def convert_all(diam, MW, gd, eng):
    np_diam = convert_np_diam(diam)
    l = convert_chain_len(MW)
    n_chains, gdens = convert_gd(diam, gd)
    eps_LJ = convert_eng(eng)
    
    ret = {'np_diam': np_diam,
           'chain_len': l,
           'n_chains': n_chains,
           'gd': gdens,
           'eps_LJ': eps_LJ
        }
    return ret
    
    
    
if __name__ == '__main__':
    
    AB_eng = 30 #kJ/mol
    gd_exp = 0.8 #0.1 to 1 
    diam_exp = 12 #8 to 35
    MW_exp = 3500 #3.5kDa & 10kDa = 26 & 73

    for gd_exp in [0.8]:
        results = convert_all(diam_exp, MW_exp, gd_exp, AB_eng)
        print(results)
    
        
    