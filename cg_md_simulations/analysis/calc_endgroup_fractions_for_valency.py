#%%
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 10:15:41 2025

@author: skronen
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
import networkx as nx
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from collections import defaultdict
import warnings
import sys
from scipy.interpolate import griddata, make_splrep
import json
import pandas as pd

"""
Some visualiazation functions are included, but the bulk of the
analysis for the paper is that included in run_analysis.
"""

def read_data(file):
    data = []
    
    sav = False
    
    with open(file, 'r') as f:
        lines = f.read().splitlines()
        
    for line in lines:
        if 'Atoms' in line:
            sav = True
        if 'Bonds' in line:
            sav = False
        
        if sav: 
            if len(line.split())!=7: continue
            foo = [float(x) for x in line.split()]
            data.append(foo)
    return np.array(data)

def make_init_vecs(data, center_subtract = False):
    ids = data[:,0]
    types = data[:,2].astype(int)
    pos = data[:,-3:]
    
    ends1 = pos[types == 5]
    ends2 = pos[types == 10]
    
    vecs1 = ends1
    vecs2 = ends2
    
    pos_copy = pos.copy()
    if center_subtract:
        center1 = pos[types ==1]
        center2 = pos[types ==6]
        pos_copy[types == 5] = pos_copy[types ==5] - center1
        pos_copy[types == 10] = pos_copy[types ==10] - center2

    col = [0 for x in vecs1]
    col.extend([1 for x in vecs2])
    
    type_incl = [5, 10]
    vecs = {ix: pos for ix, pos, t in zip(ids, pos_copy, types) if t in type_incl}
    cols = {ix : type_incl.index(t) for ix, t in zip(ids, types) if t in type_incl}

    return vecs, cols  

def compute_bond_lifetimes(bond_list, timesteps, timestep_size=1.0, normalize=True):
    global bond_presence, bond_lifetimes
    
    def normalize_bond(bond):
        return tuple(sorted(bond)) if normalize else tuple(bond)

    bond_presence = defaultdict(list)  # maps bond id to list of timesteps bond is present

    for t_idx, bonds in enumerate(bond_list):
        timestep = timesteps[t_idx]
        for bond in bonds:
            norm_bond = normalize_bond(bond)
            bond_presence[norm_bond].append(timestep)

    bond_lifetimes = defaultdict(list)  #maps bond id to list of lifetimes

    for bond, times in bond_presence.items():
        times = sorted(times)
        start = times[0]
        for i in range(1, len(times)):
            if times[i] != times[i-1] + 1:
                end = times[i-1]
                duration = (end - start + 1) * timestep_size
                bond_lifetimes[bond].append(duration)
                start = times[i]
        duration = (times[-1] - start + 1) * timestep_size
        bond_lifetimes[bond].append(duration)

    return dict(bond_lifetimes)

def generalized_gaussian(x, A, x0, sigma, n):
    return A * np.exp(-np.abs((x - x0) / sigma)**n)


def fill_nobonds(bonds, ts, tsgap, maxts):
    global ts_dict
    ts_dict = {t: b for t, b in zip(ts, bonds)}
    newts = np.arange(0, maxts + tsgap, tsgap).astype(int)
    newbonds = [ts_dict.get(tsval, []) for tsval in newts]
    return newbonds, newts

def read_bond_data(file, tsgap = 1000, maxts = 5000000):
    with open(file, 'r') as f:
            lines = f.read().splitlines()
            
    ts = []
    nums = []
    bonds = []
    bond_temp = []
    temp_ts = 0
    for i,line in enumerate(lines):
        if "TIMESTEP" in line:
            if len(bond_temp) > 0:
                bonds.append(bond_temp)
                ts.append(temp_ts)
            temp_ts = int(lines[i+1])
            bond_temp = []
        elif len(line.split()) == 1:
            continue
        else:
            spl = line.split()
            dist = float(spl[0])
            id1 = int(spl[1])
            id2 = int(spl[2])
            type1= int(spl[3])
            type2 = int(spl[4])
            app = [id1, id2]
            app.sort()
            if type1 != type2 and dist < 2.6*1.5:
                bond_temp.append(app)
    
    if max(ts) > maxts: warnings.warn('max ts is less than actual max timestep')        
    
    if len(bond_temp) > 0:
        bonds.append(bond_temp)
        ts.append(temp_ts)
        
    nums = [len(x) for x in bonds]
    
    return bonds, ts, nums

def get_bondtimes(bonds, ts, tsstep = 10):
    ret = compute_bond_lifetimes(bonds, ts)
    alltimes1 = list(ret.values())
    alltimes = np.array([x for xs in alltimes1 for x in xs])
    alltimes = alltimes * tsstep
    return alltimes

def plot_bt_distribution(nums, alltimes, ts, tsstep):
    minbin = tsstep/2
    maxbin = max(alltimes)+tsstep/2
    bins = np.arange(minbin, maxbin, tsstep)
    h,bins = np.histogram(alltimes, bins = bins, density = True)
    space = bins[1]-bins[0]
    bin_mids = [b+space/2 for b in bins[:-1]]
    plt.figure()
    plt.plot(bin_mids, h, color = 'k', label = 'raw data')
    
    p0 = [1.0, np.mean(alltimes), np.std(alltimes), 2.0] 
    bounds = ([0, -np.inf, 0, 0.5], [np.inf, np.inf, np.inf, 10]) 
    popt, _ = curve_fit(generalized_gaussian, bin_mids, h, p0=p0, bounds=bounds)

    xfit =np.linspace(minbin, maxbin , 1000)
    yfit = generalized_gaussian(xfit, *popt)
    plt.plot(xfit, yfit, colo30r = 'r', ls = '--', label = 'fit')
    plt.legend()
    plt.xlabel('bond time, timesteps')
    plt.ylabel('Probability Density')
    
    meantime = popt[1]
    meantime = np.mean(alltimes)

    avg_num = np.mean(nums)/ts[1]
    
    return meantime, avg_num

def plot_num_bonds(ts, nums):
    plt.figure()
    plt.plot(ts, nums)
    plt.xlabel('Timestep')
    plt.ylabel('Num HB')

def plot_3d_graph(edges, weights, colids, xyz_coords):
    G = nx.Graph()
    for (u, v), w in zip(edges, weights):
        G.add_edge(u, v, weight=w)

    norm = mcolors.Normalize(vmin=min(weights), vmax=max(weights))
    cmap = cm.get_cmap('coolwarm')
    edge_colors = [cmap(norm(w)) for w in weights]

    node_colors = ['red', 'blue']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    alphs = weights/max(weights)
    for (u, v), color, alpha in zip(edges, edge_colors, alphs):
        x_vals = [xyz_coords[u][0], xyz_coords[v][0]]
        y_vals = [xyz_coords[u][1], xyz_coords[v][1]]
        z_vals = [xyz_coords[u][2], xyz_coords[v][2]]
        ax.plot(x_vals, y_vals, z_vals, color=color, linewidth=2, alpha = alpha)

    for idx, (node, (x, y, z)) in enumerate(xyz_coords.items()):
        ax.scatter(x, y, z, color=node_colors[colids[node]], s=10)

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(weights)
    cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.1)
    cbar.set_label("Edge weight")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.tight_layout()
    plt.show()


def plot_proba_heatmap(xyz, probabilities):
    xyz /= np.linalg.norm(xyz, axis=1, keepdims=True)

    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    r = np.linalg.norm(xyz, axis=1)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)

    mask1 = (phi >= -np.pi/2) & (phi <= np.pi/2) #particle 1
    mask2 = ~mask1 #particle 2
    
    def wrap_phi(phi):
        return ((phi + np.pi/2) % (2 * np.pi)) - 3*np.pi/2

    def plot_group(ax, theta, phi, probs, title, wrap=False):
        if wrap:
            phi = wrap_phi(phi)

        phi_range = (-np.pi/2, np.pi/2)
        theta_lin = np.linspace(0, np.pi, 200)
        phi_lin = np.linspace(*phi_range, 200)
        phi_grid, theta_grid = np.meshgrid(phi_lin, theta_lin)

        points = np.column_stack((theta, phi))
        grid_vals = griddata(points, probs, (theta_grid, phi_grid), method='linear')
        grid_vals = np.nan_to_num(grid_vals, nan=0.0)

        im = ax.pcolormesh(phi_grid, theta_grid, grid_vals, shading='auto', cmap='inferno', vmin = 0, vmax = 1)
        ax.scatter(phi, theta, c=probs, cmap='inferno', edgecolor='white', s=20, vmin = 0, vmax = 1)
        ax.set_title(title)
        ax.set_xlabel('phi (wrapped)' if wrap else 'phi')
        ax.set_ylabel('theta')
        return im

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), constrained_layout=True)

    im1 = plot_group(axes[0], theta[mask1], phi[mask1], probabilities[mask1], 'Particle 1')

    im2 = plot_group(axes[1], theta[mask2], phi[mask2], probabilities[mask2], 'Particle 2', wrap=True)

    fig.colorbar(im1, ax=axes, label='Probability', location='right')
    plt.suptitle('Spherical Heatmaps (Split by φ range)')
    plt.show()

def gauss(x, a, b, c, d):
    return a * np.exp(-b*x**d/c)
    
def plot_radial_probs_split(pos_cs, weights, ax = None, plot = False, fit = True):
    global popt, pts
    
    if not ax and plot:
        f, ax = plt.subplots()
        
    pts = pos_cs.copy()
    
    mask = pos_cs[:,0] > 0
    split_pts = [pts[mask], pts[~mask]]
    split_wts = [weights[mask], weights[~mask]]
    
    all_angles = np.arccos(np.abs(pts[:,0]))
    all_weights = weights.copy()
    colors = ['r', 'b']
    labs = ['Particle 1', 'Particle 2']
    
    all_binmids = []
    all_smooths = []
    all_popt = []
    all_radpts = []
    for i, (pts, weights) in enumerate(zip(split_pts, split_wts, strict = True)):
        pts = pts/np.linalg.norm(pts, axis = 1, keepdims = True)
        
        dots = np.abs(pts[:,0]) #want to dot with [1,0,0] if x > 0 and [-1, 0, 0] if x < 0. This is just abs(x)
        
        angles = np.arccos(dots)   
        
        nbins = 15
        bins = np.linspace(0, np.pi/2, nbins)
        space = bins[1] - bins[0]
        bin_mids = [b + space/2 for b in bins[:-1]]
        if plot: 
            ax.scatter(angles, weights, s = 2, alpha = 0.3, color = colors[i])
        
        if fit:
            bounds = [[0, 0, 0, 0], [1, np.inf, np.inf, np.inf]]
            popt, _ = curve_fit(gauss, angles, weights, bounds = bounds, maxfev = 100000)
            bin_mids = np.linspace(0, np.pi/2, 500)
            smooth = gauss(bin_mids, *popt)
            
        elif fit:    
            #spline fitting didn't really work bc of negative values, and not fitting tail well
            sort = np.argsort(angles)
            sortang = angles[sort]
            sortweights = np.log10(weights[sort]+1e-12)
            spline = make_splrep(sortang, sortweights, s = 1)
            bin_mids = np.linspace(0, np.pi/2, 500)
            smooth = 10**spline(bin_mids)
        else:
            smooth = []
            for i, b in enumerate(bin_mids):
                mask = np.logical_and(angles > bins[i], angles < bins[i+1])
                smooth.append(np.mean(weights[mask]))
        
        
        if plot:
            ax.plot(bin_mids, smooth, color = colors[i], label = labs[i])
            ax.legend()
            ax.set_xlabel('Angle from Interparticle Axis, radians')
            ax.set_ylabel('Bond Probability')
        
        all_binmids.append(bin_mids)
        all_smooths.append(smooth)
        all_popt.append(popt)
        all_radpts.append(np.array([angles, weights]))
        
    radpts = all_radpts
    popt = all_popt
    smooth = all_smooths
    bin_mids = all_binmids
    print(popt)
    return (bin_mids, smooth), radpts   
    

def plot_radial_probs(pos_cs, weights, ax=None, plot = False, fit = True, split = False):
    
    if split:
        (bin_mids, smooth), radpts  = plot_radial_probs_split(pos_cs, weights, ax = ax, plot = plot, fit = fit)
        return (bin_mids, smooth), radpts   
    
    global popt, pts
    if not ax and plot:
        f, ax = plt.subplots()
    
    pts = pos_cs.copy()
    pts = pts/np.linalg.norm(pts, axis = 1, keepdims = True)
    
    dots = np.abs(pts[:,0]) #want to dot with [1,0,0] if x > 0 and [-1, 0, 0] if x < 0. This is just abs(x)
    
    angles = np.arccos(dots)   
    
    nbins = 15
    bins = np.linspace(0, np.pi/2, nbins)
    space = bins[1] - bins[0]
    bin_mids = [b + space/2 for b in bins[:-1]]
    if plot: 
        ax.scatter(angles, weights, s = 2, alpha = 0.3)
    
    
    
    if fit:
        bounds = [[0, 0, 0, 0], [1, np.inf, np.inf, np.inf]]
        popt, _ = curve_fit(gauss, angles, weights, bounds = bounds, maxfev = 100000)
        bin_mids = np.linspace(0, np.pi/2, 500)
        smooth = gauss(bin_mids, *popt)
        
    elif fit:    
        #spline fitting didn't really work bc of negative values, and not fitting tail well
        sort = np.argsort(angles)
        sortang = angles[sort]
        sortweights = np.log10(weights[sort]+1e-12)
        spline = make_splrep(sortang, sortweights, s = 1)
        bin_mids = np.linspace(0, np.pi/2, 500)
        smooth = 10**spline(bin_mids)
    else:
        smooth = []
        for i, b in enumerate(bin_mids):
            mask = np.logical_and(angles > bins[i], angles < bins[i+1])
            smooth.append(np.mean(weights[mask]))
    
    
    if plot:
        ax.scatter(angles, weights, s = 2, alpha = 0.2)
    
    radpts = np.vstack((angles, weights)).T
    return (bin_mids, smooth), radpts

class obj:
    def __init__(self, dname, curve, radpts):
        self.dname= dname
        self.curve = curve
        self.radpts = radpts
        self.get_params()
        
    def get_params(self):
        file = os.path.join(self.dname, 'params.json')
        with open(file, 'r') as f:
            params = json.load(f)
        self.params = params
        self.d = params['np_diam1']
        self.l = params['len_graft1']
        self.testname = self.dname.split('/')[2]
        
        
def get_proba(dname):
    file = os.path.join(dname, 'pairs.txt')
    bonds, ts, nums = read_bond_data(file)
    bonds2 = [b for b in bonds if len(b) > 0]
    allb2 = np.array([b for bs in bonds2 for b in bs])
    uniq_ids, id_cts = np.unique(allb2, return_counts = True)
    bond_freq = id_cts/len(bonds)
    return np.mean(bond_freq)


def run_analysis(dnames):
    global vecs, vecs_cs, cols, pos, pos_cs, bonds, ts, nums, weights
    foobar = []
    mws = []
    
    objs = []
    
    for dname in dnames:
        print(dname)
        file = os.path.join(dname, 'pairs.txt')
        bonds, ts, nums = read_bond_data(file)
        ts2 = np.array(range(len(ts)))
        alltimes = get_bondtimes(bonds, ts2, tsstep = ts[1])
        
        #meantime, avg_num = plot_bt_distribution(nums, alltimes, ts2, ts[1])

        #print(f'{meantime=}', f'{avg_num=}', sep = '\n')

        #plot_num_bonds(ts, nums)
        
        numbonds = [len(b) for b in bonds]
        mean_numbonds = np.mean(numbonds[300:])
    
        bonds2 = [b for b in bonds if len(b) > 0]
        allb2 = np.array([b for bs in bonds2 for b in bs])
        
        uniq, cts = np.unique(allb2, axis = 0, return_counts = True)
        
        datafile = os.path.join(dname, 'data.data')
        data = read_data(datafile)
        vecs, cols = make_init_vecs(data)
        vecs_cs, cols_cs = make_init_vecs(data, center_subtract = True)
        
        molid_map = {}
        with open(datafile) as f:
            for line in f:
                sp = line.split()
                if len(sp)==6 and sp[0].isdigit(): molid_map[int(sp[0])] = int(sp[1])
                
        # plot_3d_graph(uniq, cts, cols, vecs)
    
        uniq_ids, id_cts = np.unique(allb2, return_counts = True)
        
        pos = np.array([vecs[ix] for ix in uniq_ids])

        bond_freq = id_cts/len(bonds)
        
        weights = id_cts/len(bonds)
        print('meanprob:', np.mean(weights))
        mws.append(np.max(weights))
        print(sum(weights)/2)
        
        if plot_3d:
            norm = mcolors.Normalize(vmin=min(weights), vmax=max(weights))
            cmap = cm.get_cmap('coolwarm')
            ax = plt.figure().add_subplot(projection ='3d')
            ax.scatter(*pos.T, c = weights, cmap = cmap)
            ax.set_aspect('equal')
        
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array(weights)
            cbar = plt.colorbar(sm, ax=ax, fraction=0.02, pad=0.1)
            cbar.set_label("Fraction of Time Bonded")
            plt.tight_layout()
        
        pos_cs = np.array([vecs_cs[ix] for ix in uniq_ids])
        
        allpos_cs = np.array([x for x in vecs_cs.values()])
        mask = np.array([q in pos_cs for q in allpos_cs])
        addpos = allpos_cs[~mask]
        addwt = np.array([0]*len(addpos))
        
        pos_cs = np.concatenate((pos_cs, addpos), axis = 0)
        weights = np.append(weights, addwt)
                
        if plot_heatmaps:
            plot_proba_heatmap(pos_cs, weights)
        
        nums = np.array([len(x) for x in bonds])
        
        uniq1, cts1 = np.unique(nums, return_counts = True)
        
        if plot_numbond_dist:
            plt.figure()
            plt.plot(uniq1, cts1)
            plt.xlabel('number of bonds at a given time')
            plt.ylabel('frequency')
            
            
        foobar.append([uniq1, cts1])
        
        smooth_curve, radpt = plot_radial_probs(pos_cs, weights,ax = None, plot = False, split = SPLIT)
        
        o = obj(dname, smooth_curve, radpt)
        o.popt = popt
        objs.append(o)
        print('frac_bonds: ', mean_numbonds/len(weights)*2)
    
    return objs

    
def final_series(sig = 0.8):
    global allfracs, allallfracs, objs, o, pts
    diels = [10, 15, 22.5]
    labs = [
            '25.5nm_3.5kDa_22.5mM',
            '16.8nm_3.5kDa_9.0mM',
            '11.0nm_3.5kDa_18.8mM',
            ]
    colors = ['b', 'purple', 'r']
    markers = ['*', 'o', '^', 's']
    
    simdirs = ['final_series_gd_compare',
               f'final_series_gd{sig:0.1f}_t2',
               f'final_series_gd{sig:0.1f}_t3']
    

    allallfracs = [] 

    df_fracs = []
    df_diams = []
    df_diels = []

    for sd in simdirs:
        allfracs = []

        for diel in diels:
            nacl1 = f'sims/{sd}/sig_{sig:0.2f}_d_11.1_l_26_diel_DREP_conc_18.8/00000'
            midd1 = f'sims/{sd}/sig_{sig:0.2f}_d_16.8_l_26_diel_DREP_conc_9.0/00000'
            cscl1 = f'sims/{sd}/sig_{sig:0.2f}_d_25.5_l_26_diel_DREP_conc_22.5/00000'
            
            dnames= [test.replace('DREP', f'{diel:0.2f}') for test in [cscl1, midd1, nacl1]]
            
            objs = run_analysis(dnames)
        
            
            def invgaus(y, a, b, c, d):
                return (-c/b*np.log(y/a))**(1/d)
                
            legvar = 'testname'
            cutoff = 0.05
            cutoff_type = 'nvar'       
         
            foo = []
            for o in objs:
                popt1 = o.popt
                fracs = []
                for popt in o.popt:
                    ang = invgaus(cutoff, *popt)
                    frac = (1-np.cos(ang))/2
                    o.frac = frac
            
                    print(frac)
                    fracs.append(frac)
                foo.append(np.mean(fracs))
                df_diels.append(diel)
                df_fracs.append(np.mean(fracs))
                df_diams.append(o.d)
            allfracs.append(foo)
        allallfracs.append(allfracs)
    allallfracs = np.array(allallfracs)
    
    arr = np.vstack((df_diams, df_diels, df_fracs)).T
    df = pd.DataFrame(arr, columns = ['diam', 'diel', 'frac'])

    df.to_csv(f'../results/diamseries_{sig}.csv')

    means = np.mean(allallfracs, axis = 0)
    stds = np.std(allallfracs, axis = 0)
    
    plt.figure()
    for i, (m,s) in enumerate(zip(means.T, stds.T)):
        plt.errorbar(diels, m, s, color = colors[i], marker = markers[i], label = labs[i])
        
    plt.legend()
    plt.xlabel('Dielectric Constant')
    plt.ylabel('% of Surface Area Forming Bonds')
    plt.title('cutoff at 0.5*maximum')

def final_series_asym(sig = 0.8):
    global allfracs, allallfracs, objs, o, pts, fracs
    diels = [15]
    
    simdirs = ['asym/asym_caf2_1_gd0.8/t1',
               'asym/asym_caf2_1_gd0.8/t2',
               'asym/asym_caf2_1_gd0.8/t3']
    
    if sig == 0.5:
        simdirs = simdirs[:1]

    allallfracs = [] 

    df_fracs1 = []
    df_fracs2 = []
    df_diams = []
    df_diels = []

    for sd in simdirs:
        allfracs = []

        for diel in diels:
            caf2_1 = f'sims/{sd}/sig_{sig:0.2f}_d_12.6_l_26_diel_DREP_conc_2.3/00000'
            caf2_2 = f'sims/{sd.replace("asym_caf2_1_gd0.8", "asym_caf2_2_gd0.8")}/sig_{sig:0.2f}_d_25.5_l_26_diel_DREP_conc_7.5/00000'
            uh3_1 = f'sims/{sd.replace("asym_caf2_1_gd0.8", "asym_uh3_gd0.8")}/sig_{sig:0.2f}_d_16.8_l_73_diel_DREP_conc_1.5/00000'
            
            #diel = 22.5
            dnames= [test.replace('DREP', f'{diel:0.2f}') for test in [caf2_1, caf2_2, uh3_1]]
            
            objs = run_analysis(dnames)
        
            
            def invgaus(y, a, b, c, d):
                return (-c/b*np.log(y/a))**(1/d)
                
            legvar = 'testname'
            cutoff = 0.05
            cutoff_type = 'nvar'        
         
            foo = []
            for o in objs:
                popt1 = o.popt
                fracs = []
                for popt in o.popt:

                    ang = invgaus(cutoff, *popt)
                    frac = (1-np.cos(ang))/2
                    o.frac = frac
            
                    print(frac)
                    fracs.append(frac)
                foo.append(fracs)
                df_diels.append(diel)
                df_fracs1.append(fracs[0])
                df_fracs2.append(fracs[1])
                df_diams.append(o.d)
            allfracs.append(foo)
        allallfracs.append(allfracs)
    allallfracs = np.array(allallfracs)
    
    arr = np.vstack((df_diams, df_diels, df_fracs1, df_fracs2)).T
    df = pd.DataFrame(arr, columns = ['diam', 'diel', 'frac1', 'frac2'])

    df.to_csv(f'../results/diamseries_asym_{sig}.csv')

if __name__ == "__main__":
    plot_heatmaps = False
    plot_numbond_dist = False
    plot_3d = False
    plot_rad = False
    SPLIT = True

    final_series(0.5)
    final_series(0.8)
    final_series_asym(sig = 0.8)

    sys.exit()
