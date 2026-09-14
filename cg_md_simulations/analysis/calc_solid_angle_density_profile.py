#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 21 11:54:35 2025

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from scipy.stats import gaussian_kde
import pandas as pd
import seaborn as sns

from Dump_Reader import read_dump


def get_pos_types(file):
    dump = read_dump(file, ret_dict=True, skips_ixs=True)
    xpos = dump['xpos']
    ypos = dump['ypos']
    zpos = dump['zpos']
    pos = np.concatenate((xpos[:,:, np.newaxis], ypos[:,:, np.newaxis], zpos[:,:, np.newaxis]), axis = -1)
    types = dump['types']
    mols = dump['molecs']
    return pos, types, mols

def subtract_centers(pos, types, mols):
    mask = types > 5
    ret= []
    for m in [mask, ~mask]:
        p = pos[:,m,:]
        cmask = np.array([t in [1, 6] for t in types[m]])
        center = p[:,cmask]
        p =p - center
        p = p[:,~cmask] #remove center bead
        ret.append(p)
    return ret

def sa_cap(R, h):
    return 2*np.pi*R*h
    
def get_slice_areas(bins):
    caphts = 1*(1-np.cos(bins))
    slice_areas = []
    for i in range(len(caphts)-1):
        slicearea = sa_cap(1, caphts[i+1]) - sa_cap(1, caphts[i])
        slice_areas.append(slicearea)
    slice_areas = np.array(slice_areas)
    return slice_areas

class densprof:
    def __init__(self, dname, nbins = 100):
        self.dname = dname
        self.nbins = nbins
        self.efile = os.path.join(dname, 'equil_ends.lammpstrj')
        self.pfile = os.path.join(dname, 'ends.lammpstrj')
        self.get_dprofs()
    
    def get_dprofs(self):
        for i, file in enumerate([self.efile, self.pfile]):
            pos, types, mols = get_pos_types(file)
            pos_cs = subtract_centers(pos, types, mols)
            pos_cs = np.concatenate(pos_cs, axis = 1)
            norms = np.linalg.norm(pos_cs,axis = -1, keepdims=True)
            pos_cs = pos_cs/norms
            rad = np.mean(norms) / 10 #average radius in nm
            
            solid_angle = np.arccos(np.abs(pos_cs[:,:,0]))
            sas = solid_angle.ravel()

            bins = np.linspace(0, np.pi/2, self.nbins)
            h, bins = np.histogram(sas, bins = bins, density = False)
            space = bins[1] - bins[0]
            bin_mids = np.array([b + space/2 for b in bins[:-1]])
            
            slice_areas = get_slice_areas(bins)
            hnorm = h/slice_areas/len(pos) #units have been normalized to R =1 
            hnorm = hnorm /rad**2
            
            cumulative_counts = np.cumsum(h)

            cum_areas = np.cumsum(slice_areas)
            cum_dens = cumulative_counts/cum_areas/len(pos)/rad**2
            mask = np.logical_and(bin_mids < np.pi/4, bin_mids > 0.2) #avoid noise at high and low angles
            maxcum = np.max(cum_dens[mask])
            maxcum_r = bin_mids[np.argmax(cum_dens[bin_mids < np.pi/4])]
            if i == 0:
                self.equil_dens = np.vstack((bin_mids, hnorm))
                self.equil_cum = np.vstack((bin_mids, cum_dens))
                self.equil_maxcum = maxcum

            else:
                self.prod_dens = np.vstack((bin_mids, hnorm))
                self.prod_cum = np.vstack((bin_mids, cum_dens))
                self.prod_maxcum = maxcum

        self.ratio = self.prod_maxcum/self.equil_maxcum
        
        sa = self.prod_cum[0][np.argmax(self.prod_cum[1])]
        self.sa_perc = sa/2/np.pi
        
    def plot(self):
        plt.figure(figsize = (5, 5), dpi = 400) 
        ix2e = np.where(self.equil_dens[0] > 1.5)[0][0]
        ix2p = np.where(self.prod_dens[0] > 1.5)[0][0]

        plt.plot(*self.equil_dens[:,:ix2e], color = 'k', ls = ':', alpha = 0.6)
        plt.plot(*self.prod_dens[:,:ix2p], color = 'r', ls = ':', alpha = 0.6)
        
        cum_dens = self.equil_cum[1]
        ix = np.where(cum_dens> 0.5*np.max(cum_dens))[0][0]
        ix2 = np.where(self.equil_cum[0] > 1.5)[0][0]
        plt.plot(*self.equil_cum[:,ix:ix2], color = 'k', label = 'Isolated')
        
        cum_dens = self.prod_cum[1]
        ix = np.where(cum_dens> 0.5*np.max(cum_dens))[0][0]
        ix2 = np.where(self.prod_cum[0] > 1.5)[0][0]
        plt.plot(*self.prod_cum[:,ix:ix2], color = 'r', label = 'Contact')
        
        plt.legend()
        plt.xlabel('Solid Angle, radians')
        plt.ylabel(r'Cumulative Endgroup Density, $nm^{-2}$')
        plt.title(self.dname.split('/')[-2])
        if self.dname.split('/')[-2] == 'sig_0.80_d_25.5_l_26_diel_15.00_conc_22.5':
            plt.savefig('../paper_plots/example_densplot.png')
    
def get_ang_dens_prof(dname, plot):
    obj = densprof(dname)
    if plot:
        obj.plot()
    return obj
    
def get_dirs(folder):
    dnames = os.listdir(folder)
    dnames = [d for d in dnames if d.startswith('sig_')]
    dnames = [os.path.join(folder, d, '00000') for d in dnames]
    dnames.sort()
    return dnames
    
if __name__ == '__main__':
    folder = '/home/skronen/Documents/muri_am/mdsims/charged_endgroups_debye/sims/final_series_gd_compare'
    dnames = get_dirs(folder)
    f2 = '/home/skronen/Documents/muri_am/mdsims/charged_endgroups_debye/sims/final_series_gd0.8_t2'
    f3 = '/home/skronen/Documents/muri_am/mdsims/charged_endgroups_debye/sims/final_series_gd0.8_t3'
    dnames.extend(get_dirs(f2))
    dnames.extend(get_dirs(f3))
    
    sigs = []
    ds = []
    ls = []
    diels = []
    concs = []
    ratios = []
    equils = []
    prods = []
    
    foos = []
    ct= 0
    for i,dname in enumerate(dnames[:]):
        if not 'diel_15' in dname: continue
        if not 'sig_0.80' in dname: continue
        if not 'd_25.5_l_26' in dname : continue
        spl = dname.split('/')[-2].split('_')
        sig = float(spl[1])
        #if sig !=0.8: continue
        d = float(spl[3])
        l = float(spl[5])
        diel = float(spl[7])
        conc = float(spl[9])
        
        sigs.append(sig)
        ds.append(d)
        ls.append(l)
        diels.append(diel)
        concs.append(conc)
        if ct == 8:
            plot =True
        else:
            plot = True
        foo = get_ang_dens_prof(dname, plot = plot) 
        if ct > 0:    
            sys.exit()
        ratio = foo.ratio
        prodmax = foo.prod_maxcum
        equilmax = foo.equil_maxcum
        
        foos.append(foo)
        print(f'{ratio=}')
        ratios.append(ratio)
        equils.append(equilmax)
        prods.append(prodmax)
        ct+=1
        
    percs = [a.sa_perc for a in foos]
    arr = np.vstack((sigs, ds, ls, diels, concs, ratios, equils, prods, percs)).T
    colnames = ['sig', 'd', 'l', 'diel', 'conc', 'ratio', 'equil', 'prod', 'perc']
    df = pd.DataFrame(arr, columns = colnames)
    df.to_csv('../results/density_ratios.csv')
