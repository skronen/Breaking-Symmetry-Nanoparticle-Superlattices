#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:26:25 2024

@author: skronen
"""

import os
import subprocess
import shutil
import sys

from debye_length_calc import convert_conc_to_debye


Ds = [11.1, 16.8, 25.5] #nanoparticle core diameters (nm)
sigs = [0.8] #grafting densities (nm^-2)
Ls = [26] #graft length (number of beads)
diels = [10, 15, 22.5] #dielectric constant
concs = [18.8, 9, 22.5] #NaCl concentration (mM)

series_name = 'sims/final_series_gd0.8_t1'

debyes = convert_conc_to_debye(concs)
print(debyes)


if os.path.isdir(series_name):
    shutil.rmtree(series_name)
    
def get_indexes(folder):
    ifile = os.path.join(folder, 'index.txt')
    with open(ifile, 'r') as f:
        ret = f.read().splitlines()
    return ret
    
def read_all_series_indexes():
    folders = os.listdir(series_name)
    all_index = []
    for folder in folders:
        index = get_indexes(os.path.join(series_name, folder))
        index = [os.path.join(folder, i) for i in index]
        all_index.extend(index)
    
    lines = [i + '\n' for i in all_index]
    ifile = os.path.join(series_name, 'index.txt')
    
    with open(ifile, 'w') as f:
        f.writelines(lines)
    return len(lines)


def add_num_arraytasks(array_script, N):
    if N > 90: 
        N = str(N) + '%90'
    else:
        N = str(N)
    with open(array_script, 'r') as f:
        lines = f.read().splitlines()
    
    with open(array_script, 'w') as f:
        for line in lines:
            if line == '#SBATCH --array=1-N_ARRAY_TASKS_REPLACE':
                line = line.replace('N_ARRAY_TASKS_REPLACE', N)
            f.write(line + '\n')
            

for D, debye, conc in zip(Ds, debyes, concs):
    for sig in sigs:
        for L in Ls:
            for diel in diels:
                print(D, sig, L, diel)
                callstring = ['python3', 'make_series.py', series_name, str(D), str(sig), str(L), str(diel), str(debye), str(conc)]
                subprocess.call(callstring)
                
                
N = read_all_series_indexes()
print(N)
array_script = 'array_sub.qs'
shutil.copy(array_script, series_name)
add_num_arraytasks(os.path.join(series_name, array_script), N)
    
sys.exit()




