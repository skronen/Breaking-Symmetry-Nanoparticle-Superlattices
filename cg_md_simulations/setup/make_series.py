#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 15:01:40 2024

@author: skronen
"""

import sys
import os
import subprocess
import shutil
import numpy as np
import json
from debye_length_calc import convert_conc_to_debye


D = float(sys.argv[2])
sig = float(sys.argv[3])
L = int(sys.argv[4])
diel = float(sys.argv[5])
debye = float(sys.argv[6])
conc= float(sys.argv[7])
sname = sys.argv[1]
D1 = D2 = D
L1 = L2 = L

#For asymmetric systems, need to define separate NCT designs here, as below.

# D1 = 12.6
# D2 = 12.6
# L1 = 26
# L2 = 73
# sig = 0.8
# conc = 2.3
# debye = convert_conc_to_debye(conc)
# diel = 15
# sname = 'sims/asym_caf2_1_gd0.8/t3'

# D1 = 25.5
# D2 = 12.6
# L1 = 26
# L2 = 26
# sig = 0.8
# conc = 7.5 
# debye = convert_conc_to_debye(conc)
# diel = 15
# sname = 'sims/asym_caf2_2_gd0.8/t3'

# D1 = 16.8
# D2 = 12.6
# L1 = 73
# L2 = 26
# sig = 0.8
# conc = 1.5 
# debye = convert_conc_to_debye(conc)
# diel = 15
# sname = 'sims/asym_uh3_gd0.8/t3'


seriesname = f'sig_{sig:0.2f}_d_{D1:02}_l_{L1:02}_diel_{diel:0.2f}_conc_{conc}'
seriesname = os.path.join(sname, seriesname)


def write_json(vals, fname = 'params.json'):
    global json_dict
    
    names = ['num_np1', 'num_np2', 'num_graft1', 'num_graft2',
        'len_graft1', 'len_graft2', 'frac_sticky1', 'frac_sticky2',
        'np_diam1', 'np_diam2', 'graft_diam', 'graft_eps', 'opp_sticky_eps',
        'same_sticky_eps', 'sticky_embed_distance', 'sticky_diam', 'graft_stiffness_k (NOTUSED)', 'dielectric', 'debye_length (angstroms)',
        'graft_density1', 'graft_density2',]
     
    json_dict = {name:val for name, val in zip(names, vals)}
    with open(fname, 'w') as f:
        json.dump(json_dict, f, indent = 4)

def add_num_arraytasks(array_script):
    with open(array_script, 'r') as f:
        lines = f.read().splitlines()
    
    with open(array_script, 'w') as f:
        for line in lines:
            if line == '#SBATCH --array=1-N_ARRAY_TASKS_REPLACE':
                line = line.replace('N_ARRAY_TASKS_REPLACE', str(len(dirlist)))
            
            f.write(line + '\n')

def write_radfile(dirname):
    file = os.path.join(dirname, 'radius.txt')
    with open(file, 'w') as f:
        f.write(str((arglist[8] + arglist[9])/2))

def convert_gd(diam, gd, sigma= 1):
    SA = np.pi * diam**2 / 100 #convert to angstroms
    n_chains = int(np.round(gd * SA, 0))
    gdens = gd * sigma**2
    return n_chains, gdens
        
if __name__ == '__main__':
    trials = 1
    n_np1 = [1]
    n_np2 = [1]
    gd1 = [sig]
    gd2 = [sig]
    len_graft1 = [L1]
    len_graft2 = [L2]
    frac1 = [1] 
    frac2 = [1] 
    np_diam1 = [D1*10]
    np_diam2 = [D2*10]
    graft_diam = [7.6]
    graft_eps = [0]
    op_sticky_eps = [13]
    same_sticky_eps = [1]
    sticky_embed_distance = [6.3] #1 is fully embedded, 0 is not embedded
    sticky_diam = [0.3*7.6]
    graft_k = [0.9] #not actually used for these charged sims
    dielectric = [diel]
    debyes = [debye]
    
    (gd1, gd2, len_graft1, len_graft2, frac1, frac2,
        np_diam1, np_diam2, graft_diam, graft_eps, op_sticky_eps,
            sticky_embed_distance, sticky_diam, graft_k) = [i *len(same_sticky_eps) for i in (gd1, gd2, len_graft1, len_graft2, frac1, frac2,
                    np_diam1, np_diam2, graft_diam, graft_eps, op_sticky_eps,
                        sticky_embed_distance, sticky_diam, graft_k)]

    
    if trials > 1:
        (gd1, gd2, len_graft1, len_graft2, frac1, frac2,
            np_diam1, np_diam2, graft_diam, graft_eps, op_sticky_eps,
                same_sticky_eps, sticky_embed_distance, sticky_diam, graft_k, dielectric, debyes) = [i * trials for i in (gd1, gd2, len_graft1, len_graft2, frac1, frac2,
                        np_diam1, np_diam2, graft_diam, graft_eps, op_sticky_eps,
                            same_sticky_eps, sticky_embed_distance, sticky_diam, graft_k)]
                                                                                             
    dirlist = []
    for i in range(len(gd1)):
        runlist = ['python3', 'make_pgp_data_onetoone.py']
        
        
        ng_1, _= convert_gd(np_diam1[i], gd1[i], sigma = graft_diam[i])
        ng_2, _= convert_gd(np_diam2[i], gd2[i], sigma = graft_diam[i])

        arglist = [1, 
                  1, 
                  ng_1, 
                  ng_2, 
                  len_graft1[i],
                  len_graft2[i],
                  frac1[i],
                  frac2[i],
                  np_diam1[i],
                  np_diam2[i], 
                  graft_diam[i], 
                  graft_eps[i], 
                  op_sticky_eps[i],
                  same_sticky_eps[i],
                  sticky_embed_distance[i], 
                  sticky_diam[i],
                  graft_k[i],
                  dielectric[i],
                  debyes[i]]
        
        
        #CHECK THAT GRAFTS WONT OVERLAP:
        working_ratio = 1.5 #maximum ratio that I know runs alright
        sa_np1 = 4*np.pi*(np_diam1[i]/2)**2
        sa_graft1 = np.pi*(graft_diam[i]/2)**2*ng_1
        rat1 = sa_graft1/sa_np1
        if rat1>working_ratio:
            print(f'Warning sim {i}: grafts may be too dense on particle 1 ({rat1})')
        sa_np2 = 4*np.pi*(np_diam2[i]/2)**2
        sa_graft2 = np.pi*(graft_diam[i]/2)**2*ng_2
        rat2 = sa_graft2/sa_np2
        if rat2>working_ratio:
            print(f'Warning sim {i}: grafts may be too dense on particle 1 ({rat2})')
        
        arglist.append(gd1[i])
        arglist.append(gd2[i])
        
        runlist.extend([str(a) for a in arglist])
        
        dirname = f'{seriesname}/{(i):05}'
        dirlist.append(dirname.split('/')[-1])
        try:
            os.makedirs(dirname)
        except:pass
        out = subprocess.call(runlist)
        
        write_json(arglist)
        for f in ['data.data', 'pgp_equil.in', 'params.json', 'pgp_prod.in']:
            if os.path.exists(os.path.join(dirname, f)):
                os.remove(os.path.join(dirname, f))
            shutil.move(f, dirname)
            
        shutil.copy('cleanup_rees.sh', dirname)
        shutil.copy('filter_bonds.sh', dirname)

        write_radfile(dirname)
            
        for jname in ['array_sub' ]:
            shutil.copy(f'{jname}.qs', seriesname)
        
        add_num_arraytasks(os.path.join(seriesname, 'array_sub.qs'))
        
        assert len(dirlist) < 50
        
        with open(f'{seriesname}/index.txt', 'w') as f:
            for d in dirlist:
                f.write(f'{d}\n')
                
        
                        
