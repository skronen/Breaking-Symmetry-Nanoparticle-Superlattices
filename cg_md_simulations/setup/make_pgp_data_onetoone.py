#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:30:22 2024

@author: skronen
"""

import numpy as np
import math
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt
import sys

COUL_CUTOFF = 140


fene = True
endbead_size = 2.6 #angstroms

"""params to vary"""
N_pgp = [int(sys.argv[1]), int(sys.argv[2])] #number of PGPs
N_grafts = [int(sys.argv[3]), int(sys.argv[4])] # Number of Grafts
len_grafts = [int(sys.argv[5]), int(sys.argv[6])] # length of grafts

np_size = [float(sys.argv[9]), float(sys.argv[10])] #np diam size
graft_bead_size = float(sys.argv[11]) #size of graft chain beads
graft_eps = float(sys.argv[12])
ang_k_graft = float(sys.argv[17])
dielectric = float(sys.argv[18])
debye_length = float(sys.argv[19])

if not endbead_size:
    endbead_size = graft_bead_size

template_equil = 'equil.in'
template_prod = 'prod.in'

          
"""sticky params"""
sticky_fracs = [float(sys.argv[7]), float(sys.argv[8])]
if np.any(sticky_fracs):
    sticky_ends = True
else:sticky_ends = False

sticky_eps = float(sys.argv[13])
same_sticky_eps = float(sys.argv[14])
embed_distance = float(sys.argv[15]) #distance at which to embed sticky bead (fraction of graft_bead diameter)-- 0 is 1 graft bead diameter away
sig_sticky = float(sys.argv[16])

""" constant params"""

if graft_eps == 0: graft_wca = True
else:graft_wca = False
    
if sticky_ends:
    types_per_np = 5
else:
    types_per_np = 3
    
def calc_init_len_box():
    diams = np.array(np_size) + 2* np.array(len_grafts)*graft_bead_size
    distance = np.sum(diams)/2 + 5*debye_length #3 is outside cutoff of 2.5 sigma
    space_factor = 1.5
    len_box = 0.5 * (np.sum(diams/2) + distance ) * space_factor
    len_box2 = 0.5 * np.max(diams) * space_factor
    return distance, len_box, len_box2

distance, _, len_box2 = calc_init_len_box()
no_pgptypes = 2
len_box = distance + 0.02
final_box = len_box
bond_k = 30
ang_k = 50

    
    
    
def make_np(diam):
    sa_np = 4 * np.pi * (diam/2)**2
    a_bead = np.pi* (0.5*graft_bead_size) **2
    num_beads = int(sa_np/a_bead *2) #multiply by 2 to ensure coverage
    points = fibonacci_sphere(samples = num_beads)
    points *= diam/2
    points = np.append([[0,0,0]], points, axis = 0) #put center bead at front
    types = [1]
    types.extend([2 for i in range(num_beads)])
    return points, types

def make_nb_params(graft_wca = True):
    sigs = []
    nb_params = []
    
    if graft_wca:
        for i in range(len(N_pgp)):
            tb = types_per_np * i #type_base
            sig11 = 1
            sig13 = 1/2 + graft_bead_size/2
            sig33 = graft_bead_size
            sig15 = (sig11 + endbead_size)/2
            sig35 = (sig33 + endbead_size)/2
            sig55 = endbead_size

            params = [[tb+1, tb+1, 1, sig11, 2**(1/6) * sig11],
                      [tb+1, tb+2, 1, sig13, 2**(1/6) * sig13],
                      [tb+2, tb+2, 1, sig33, 2**(1/6) * sig33],
                      [tb+2, tb+3, 1, sig33, 2**(1/6) * sig33],
                      [tb+1, tb+3, 1, sig13, 2**(1/6) * sig13],
                      [tb+3, tb+3, 1, sig33, 2**(1/6) * sig33]]
            nb_params.extend(params)
            
            if sticky_ends:
                sig14 = sig_sticky/2 + 1/2
                sig34 = sig_sticky/2 + graft_bead_size/2
                sig44 = sig_sticky
                sig45 = sig_sticky/2 + endbead_size/2
                app =  [[tb+1, tb+4, 1, sig14, 2**(1/6) * sig14],
                        [tb+2, tb+4, 1, sig34, 2**(1/6) * sig34],
                        [tb+3, tb+4, 1, sig34, 2**(1/6) * sig34],
                        [tb+4, tb+4, 1, 2.3*sig44, 2.3*2**(1/6) * sig44],
                        [tb+1, tb+5, 1, sig15, 2**(1/6) * sig15], 
                        [tb+2, tb+5, 1, sig35, 2**(1/6) * sig33], 
                        [tb+3, tb+5, 1, sig35, 2**(1/6) * sig35], 
                        [tb+4, tb+5, 1, sig45, 2**(1/6) * sig34], 
                        [tb+5, tb+5, 1, sig55, 2**(1/6) * sig55, COUL_CUTOFF]
                        ]
                nb_params.extend(app)
                sigs.append([sig11, sig33, sig44])
            else:
                sigs.append([sig11, sig33])
    else:
        for i in range(len(N_pgp)):
            tb = types_per_np * i #type_base
            sig11 = 1
            sig13 = 1/2 + graft_bead_size/2
            sig33 = graft_bead_size
            sig15 = (sig11 + endbead_size)/2
            sig35 = (sig33 + endbead_size)/2
            sig55 = endbead_size
            params = [[tb+1, tb+1, 1, sig11, 2**(1/6) * sig11],
                      [tb+1, tb+2, 1, sig13, 2**(1/6) * sig13],
                      [tb+2, tb+2, 1, sig33, 2**(1/6) * sig33],
                      [tb+2, tb+3, 1, sig33, 2**(1/6) * sig33],
                      [tb+1, tb+3, 1, sig13, 2**(1/6) * sig13],
                      [tb+3, tb+3, graft_eps, sig33]]
            nb_params.extend(params)
            
            if sticky_ends:
                sig14 = sig_sticky/2 + 1/2
                sig34 = sig_sticky/2 + graft_bead_size/2
                sig44 = sig_sticky
                sig45 = sig_sticky/2 + endbead_size/2
                app =  [[tb+1, tb+4, 1, sig14, 2**(1/6) * sig14],
                        [tb+2, tb+4, 1, sig34, 2**(1/6) * sig34],
                        [tb+3, tb+4, 1, sig34, 2**(1/6) * sig34],
                        [tb+4, tb+4, 1, 2.3*sig44, 2.3*2**(1/6) * sig44],
                        [tb+1, tb+5, 1, sig15, 2**(1/6) * sig15], 
                        [tb+2, tb+5, 1, sig35, 2**(1/6) * sig35], 
                        [tb+3, tb+5, 1, sig35, 2**(1/6) * sig35], 
                        [tb+4, tb+5, 1, sig45, 2**(1/6) * sig45], 
                        [tb+5, tb+5, 1, sig55, 2**(1/6) * sig55, COUL_CUTOFF]
                        ]
                nb_params.extend(app)
                sigs.append([sig11, sig33, sig44])
            else:
                sigs.append([sig11, sig33])
            
    sigs = np.array(sigs)
    
    if sticky_ends:
        cross_params = [[1, 6, 1, sig11, 2**(1/6) * sig11],
                        [1, 7, 1, sig13, 2**(1/6) * sig13],
                        [1, 8, 1, sig13, 2**(1/6) * sig13],
                        [1, 9, 1, sig14, 2**(1/6) * sig14],
                        [2, 6, 1, sig13, 2**(1/6) * sig13],
                        [2, 7, 1, sig33, 2**(1/6) * sig33],
                        [2, 8, 1, sig33, 2**(1/6) * sig33],
                        [2, 9, 1, sig34, 2**(1/6) * sig34],
                        [3, 6, 1, sig13, 2**(1/6) * sig13],
                        [3, 7, 1, sig33, 2**(1/6) * sig33],

                        [3, 9, 1, sig34, 2**(1/6) * sig34],
                        [4, 6, 1, sig14, 2**(1/6) * sig14],
                        [4, 7, 1, sig34, 2**(1/6) * sig14],
                        [4, 8, 1, sig34, 2**(1/6) * sig34],
                        [4, 9, sticky_eps, sig44],
                        
                        [5, 6, 1, sig15, 2**(1/6) * sig15],
                        [5, 7, 1, sig35, 2**(1/6) * sig35],
                        [5, 8, 1, sig35, 2**(1/6) * sig35],
                        [5, 9, 1, sig45, 2**(1/6) * sig45],
                        [5, 10, 1, sig55, 2**(1/6) * sig55, COUL_CUTOFF],
                        
                        [1, 10, 1, sig15, 2**(1/6) * sig15],
                        [2, 10, 1, sig35, 2**(1/6) * sig35],
                        [3, 10, 1, sig35, 2**(1/6) * sig35],
                        [4, 10, 1, sig45, 2**(1/6) * sig45],
                        ]
        if graft_wca: 
            cross_params.append([3, 8, 1, sig33, 2**(1/6) * sig33])
        else:
            cross_params.append([3, 8, graft_eps, sig33])
                        
    else:
        cross_params = [[1, 4, 1, sig11, 2**(1/6) * sig11],
                        [1, 5, 1, sig15, 2**(1/6) * sig15],
                        [1, 6, 1, sig13, 2**(1/6) * sig13],
                        [2, 4, 1, sig11, 2**(1/6) * sig11],
                        [2, 5, 1, sig35, 2**(1/6) * sig35],
                        [2, 6, 1, sig13, 2**(1/6) * sig13],
                        [3, 4, 1, sig13, 2**(1/6) * sig13],
                        [3, 5, 1, sig35, 2**(1/6) * sig35]]
                        
        if graft_wca: 
            cross_params.append([3, 6, 1, sig33, 2**(1/6) * sig33])
        else:
            cross_params.append([3, 6, graft_eps, sig33])

    nb_params.extend(cross_params)
    return nb_params

def make_nb_params_equil():
    global sigs
    sigs = []
    nb_params = []

    for i in range(2):
        tb = types_per_np * i #type_base
        sig11 = 1
        sig13 = 1/2 + graft_bead_size/2
        sig33 = graft_bead_size
        sig15 = (sig11 + endbead_size)/2
        sig35 = (sig33 + endbead_size)/2
        sig55 = endbead_size
        
        params = [[tb+1, tb+1, 1, sig11, 2**(1/6) * sig11],
                  [tb+1, tb+2, 1, sig13, 2**(1/6) * sig13],
                  [tb+2, tb+2, 1, sig33, 2**(1/6) * sig33],
                  [tb+2, tb+3, 1, sig33, 2**(1/6) * sig33],
                  [tb+1, tb+3, 1, sig13, 2**(1/6) * sig13],
                  [tb+3, tb+3, 1, sig33, 2**(1/6) * sig33],
                  ]
        nb_params.extend(params)
    
    
    
        if sticky_ends:
            sig14 = sig_sticky/2 + 1/2
            sig34 = sig_sticky/2 + graft_bead_size/2
            sig44 = sig_sticky
            sig45 = sig_sticky/2 + endbead_size/2
            app =  [[tb+1, tb+4, 1, sig14, 2**(1/6) * sig14],
                    [tb+2, tb+4, 1, sig34, 2**(1/6) * sig34],
                    [tb+3, tb+4, 1, sig34, 2**(1/6) * sig34],
                    [tb+4, tb+4, 1, 2.3*sig44, 2.3*2**(1/6) * sig44],
                    [tb+1, tb+5, 1, sig15, 2**(1/6) * sig15], 
                    [tb+2, tb+5, 1, sig35, 2**(1/6) * sig35], 
                    [tb+3, tb+5, 1, sig35, 2**(1/6) * sig35], 
                    [tb+4, tb+5, 1, sig45, 2**(1/6) * sig45], 
                    [tb+5, tb+5, 1, sig55, 2**(1/6) * sig55, COUL_CUTOFF],
                    ]
            nb_params.extend(app)
            sigs.append([sig11, sig33, sig44])
        else:
            sigs.append([sig11, sig33])
            
    sigs = np.array(sigs)
    
    if sticky_ends:
        cross_params = [[1, 6, 1, sig11, 2**(1/6) * sig11],
                        [1, 7, 1, sig13, 2**(1/6) * sig13],
                        [1, 8, 1, sig13, 2**(1/6) * sig13],
                        [1, 9, 1, sig14, 2**(1/6) * sig14],
                        [2, 6, 1, sig13, 2**(1/6) * sig13],
                        [2, 7, 1, sig33, 2**(1/6) * sig33],
                        [2, 8, 1, sig33, 2**(1/6) * sig33],
                        [2, 9, 1, sig34, 2**(1/6) * sig34],
                        [3, 6, 1, sig13, 2**(1/6) * sig13],
                        [3, 7, 1, sig33, 2**(1/6) * sig33],
                        [3, 8, 1, sig33, 2**(1/6) * sig33],
                        [3, 9, 1, sig34, 2**(1/6) * sig34],
                        [4, 6, 1, sig14, 2**(1/6) * sig14],
                        [4, 7, 1, sig34, 2**(1/6) * sig34],
                        [4, 8, 1, sig34, 2**(1/6) * sig34],
                        [4, 9, 1, sig44, 2**(1/6) * sig44],
        
                        [5, 6, 1, sig15, 2**(1/6) * sig15],
                        [5, 7, 1, sig35, 2**(1/6) * sig35],
                        [5, 8, 1, sig35, 2**(1/6) * sig35],
                        [5, 9, 1, sig45, 2**(1/6) * sig45],
                        
                        [1, 10, 1, sig15, 2**(1/6) * sig15],
                        [2, 10, 1, sig35, 2**(1/6) * sig35],
                        [3, 10, 1, sig35, 2**(1/6) * sig35],
                        [4, 10, 1, sig45, 2**(1/6) * sig45],
                        
                        [5, 10, 1, sig55, 2**(1/6) * sig55, COUL_CUTOFF],

                        ]
        
    else:
        cross_params = [[1, 4, 1, sig11, 2**(1/6) * sig11],
                        [1, 5, 1, sig11, 2**(1/6) * sig11],
                        [1, 6, 1, sig13, 2**(1/6) * sig13],
                        [2, 4, 1, sig11, 2**(1/6) * sig11],
                        [2, 5, 1, sig11, 2**(1/6) * sig11],
                        [2, 6, 1, sig13, 2**(1/6) * sig13],
                        [3, 4, 1, sig13, 2**(1/6) * sig13],
                        [3, 5, 1, sig13, 2**(1/6) * sig13],
                        [3, 6, 1, sig33, 2**(1/6) * sig33]]
                       
        
    nb_params.extend(cross_params)
    
    
    return nb_params

def make_bond_params():
    b_params = []
    ct = 1
    
    b_params.append([ct, bond_k, graft_bead_size])
    ct+=1
    if sticky_ends:
        size = (graft_bead_size+endbead_size)/2
        b_params.append([ct, bond_k, size])
        ct+=1

    return np.array(b_params)
            
        
def make_ang_params():
    a_params = []
    a_params.append([1, ang_k, 180])
    a_params.append([2, ang_k_graft, 180])

    return np.array(a_params)

def generate_random_rotation():
    random_quaternion = np.random.rand(4)
    random_quaternion /= np.linalg.norm(random_quaternion)
    
    rot = Rotation.from_quat(random_quaternion)
    return rot

def fibonacci_sphere(samples=1000, direction = [0,1,0], rot = 180):
    """
    Gets direction vectors reletively evenly spaced over the surface of
    the unit sphere using the fibonnaci sphere. These are used
    as directions for the q-vectors, and are averaged over in the scattering calculation
    to improve sampling.
    
    Parameters
    ----------
    samples : int, optional
        The number of sampled vectors. With previous testing, I've found the scattering
        profile changes if you use less than 300, but above 300, it becomes fairly constant.
        The default is 1000.
    direction : 1x3 list, optional
        if samples == 1, this direction is used. The default is [0,1,0].

    Returns
    -------
    points : n_samples x 3 array 
        Array of sampled unit vectors

    """
    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
    
    if samples == 1:
        print(direction)
        point = np.array([direction])
        point = point/np.linalg.norm(point, axis = 1)
        return point
    
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))
    points = np.array(points)
    return points

def make_origins():
    origin = np.array([[-distance/2, 0, 0], [distance/2, 0, 0]])
    return origin

def make_grouplines(groups):
    lines = []
    gnames = []
    i = 1
    for i,group in enumerate(groups):
        gpname = f'gp{i+1}'
        line1 = f'group {gpname} id '
        for g in group:
            line1 = line1 + f'{g} '
        line2 = f'neigh_modify exclude molecule/intra {gpname}'
        lines.append(line1)
        lines.append(line2)
        gnames.append(gpname)
    return lines, gnames

def calc_init_pos1():
    fraction_heuristic = 1
    pos1 = -distance/2 * fraction_heuristic
    return pos1


def calc_init_pos2():
    xmin = distance/2 - np_size[0]/2 - np_size[1]/2
    return xmin

def file_replace(filename, replace_strings, replace_vals):
    with open(filename, 'r') as f:
        lines = f.read().splitlines()
    
    finlines = []
        
    for line in lines:
        for string, val in zip(replace_strings, replace_vals, strict = True):
            if string in line:
                #print(string, line)
                line = line.replace(string, val, 1)
                #print(line)
        finlines.append(line + '\n')
        
    with open(filename, 'w') as f:
        f.writelines(finlines)  
        

def do_setup():
    #GENERATE INITIAL POSITIONS
    pts = []
    bonds = []
    if sticky_ends:
        angles = []
    types = []
    mols = []
    centers = make_origins()

    xmin = centers[0][0] - 0.01
    xmax = centers[1][0] + 0.01
    
    ct = 0 #counter for number of pgp types
    id_counter = 0 #counter for index of beads as they are added
    groups = []
    for i in range(no_pgptypes): #loop over different number of pgp types
        for i_pgp in range(N_pgp[i]): #loop over number of pgps
            sticky_counter = 0
            group = []
            pgp_pts = []
            
            np_pts, np_types = make_np(np_size[i])
            np_types = [t + i * types_per_np for t in np_types]
            np_types_crop = []
            ct1 = 0
            for pt, typ in zip(np_pts, np_types, strict = True):
                if pt[0]+centers[ct][0] < xmin or pt[0]+centers[ct][0] > xmax: 
                    continue
                pgp_pts.append(pt)
                np_types_crop.append(typ)
                id_counter +=1
                if ct1 ==0:
                    center_id = id_counter
                    ct1+=1
                groups.append(id_counter)
            np_types = np_types_crop

            
            graft_vecs = fibonacci_sphere(samples = N_grafts[i]) #sample directions for grafts from fibonacci sphere
            
            if i == 1: #rotate second PGP so that they have the same number of end groups and charge is balanced
                rot = Rotation.from_euler('z', 180, degrees = True)
                graft_vecs = rot.apply(graft_vecs)
                
            
            for vec in graft_vecs: #loop over graft directions
                graft_pt = (np_size[i]/2 + graft_bead_size/2) * vec
                if graft_pt[0]+centers[ct][0] < xmin or graft_pt[0]+centers[ct][0] > xmax: 
                    continue
                np_types.append(types_per_np*(i)+3)
                pgp_pts.append(graft_pt)
                id_counter +=1
                groups.append(id_counter)
                


                graft_ids = []
                for j in range(len_grafts[i]-1): #loop over beads in the current graft
                    pgp_pts.append(graft_pt + (j+1) * graft_bead_size * vec)
                    id_counter += 1
                    graft_ids.append(id_counter)
                    
                    if not j == len_grafts[i] -2: #if not at end of chain
                        np_types.append(types_per_np*(i)+3)

                    
                if sticky_ends:
                    rand_choice = np.random.rand() < sticky_fracs[i]
                    if rand_choice:
                        np_types.append(types_per_np*(i)+5)
                        sticky_counter +=1

                    else:
                        np_types.append(types_per_np*(i)+3)
                    
                for ixg, g in enumerate(graft_ids, 1):
                    if ixg == len(graft_ids):
                        bt = 2
                    else:
                        bt = 1
                    bonds.append([g, g-1, bt]) #index1 index2 type
                    
                for j,g in enumerate(graft_ids[:]): #start at id1 to keep all angles within a chain
                    if j ==0:  #add stiffness to the center atom so that chains stay relatively perpendicular to surface, when high persistence
                        angles.append([g, g-1, center_id, 2]) #index1 index2 index3 type
                    else:
                        angles.append([g, g-1, g-2, 2]) #index1 index2 index3 type

            #handle atom types
            pgp_types = np_types
            pts.append(np.array(pgp_pts) + centers[ct])
            types.append(pgp_types)
            mols.append((i * N_pgp[0]) + i_pgp+1)
            ct+=1    
    
    #PLOT RESULTING POSTIONS
    plot = False
    if plot:
        colors = ['r','orange', 'gold', 'g', 'b', 'purple']
        
        ax = plt.figure().add_subplot(projection = '3d')
        for pt, typ in zip(pts, types):
            c = [colors[t-1] for t in typ]
            ax.scatter(pt.T[0], pt.T[1], pt.T[2], color = c)
        ax.set_aspect('equal')
        
        
    box_str = ''
    box_size = [[-len_box/2-0.49,len_box/2+0.49],[-len_box2, len_box2],[-len_box2, len_box2]]
    for dim,i in (zip(['x','y','z'],range(3))):
        for j in (range(2)):
            box_str+= (str(box_size[i][j]) + ' ')
        box_str += (dim+'lo '+dim+'hi\n')
    
    nb_params = make_nb_params_equil()
    b_params = make_bond_params()
    if sticky_ends:
        a_params = make_ang_params()
    
    noat_t = types_per_np * no_pgptypes
    nob_t = len(b_params)
    noan_t = 2
    
    num_types = no_pgptypes * types_per_np
    
    flat_types = [x for xs in types for x in xs]
    uniq, nums = np.unique(flat_types, return_counts = True)
    print(nums[uniq==5], nums[uniq==10])
    
    
    noat = np.sum([len(pt) for pt in pts] ) 

    nob = len(bonds) 
    if sticky_ends:
        noan = len(angles)
    
    #WRITE DATAFILE
    if sticky_ends:
        masses = np.array([1 if (i!=0 and i!=5) else 0.01 for i in range(num_types)])
    else:
        masses = np.array([1 if (i!=0 and i!=3) else 0.01 for i in range(num_types)])
    filename = 'data.data'

    with open(filename,'w') as F:
        F.write('Polymer-grafted-nanoparticle Data\n\n')
                    
        F.write(str(noat)+' atoms\n')
        F.write(str(nob)+' bonds\n')
        if sticky_ends:
            F.write(str(noan)+' angles\n')
        
        F.write(str(noat_t)+' atom types\n')
        F.write(str(nob_t)+' bond types\n')
        if sticky_ends:
            F.write(str(noan_t)+' angle types\n')
        
        F.write(box_str)
        
        F.write('\nMasses\n\n')
        for i, mass in enumerate(masses):
            F.write(f'{i+1} {mass}\n')
    
        F.write('\nAtoms\n\n')
        ind = 1
        for pt, typ, mol in zip(pts, types, mols, strict = True):
            for p, t in zip(pt, typ, strict = True):
                if t == 5:
                    q = 1
                elif t == 10:
                    q = -1
                else: 
                    q= 0
                F.write(f'{ind} {mol} {t} {q} {p[0]} {p[1]} {p[2]}\n')
                ind +=1
    
            
        F.write('\nBonds\n\n')
        ind = 1
        for b in bonds:
            F.write(f'{ind} {b[2]} {b[0]} {b[1]}\n')
            ind +=1
        
        ind = 1
        if sticky_ends:
            F.write('\nAngles\n\n')
            for a in angles:
                F.write(f'{ind} {a[3]} {a[0]} {a[1]} {a[2]}\n')
                ind +=1
            
    equil_file = 'pgp_equil.in'
    with open(template_equil, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    
    #WRITE INPUT EQUIL SCRIPT
    write_groups= True #ensure you only write the group lines once

    with open(equil_file,'w') as f:
        for line in lines:
            if 'box_final' in line:
                line = line.replace('50', str(final_box))
            if line=='variable rseed equal RSEEDREPLACE':
                rseed = int(np.random.rand()*10000000)
                line=f'variable rseed equal \"{rseed}\"'
            
            if line=='group g2 type 4' and sticky_ends:
                line = 'group g2 type 5'
            
            if line=='group centers type 1 4' and sticky_ends:
                line='group centers type 1 5'
            
            if fene and line == 'bond_style harmonic':
                line = 'bond_style hybrid harmonic fene'
            
            if 'DIELREPLACE' in line:
                line = line.replace('DIELREPLACE', str(dielectric))
            if 'DEBYEREPLACE' in line:
                line = line.replace('DEBYEREPLACE', str(1/debye_length)) #take inverse of debye length here
            if 'COULCUTREPLACE' in line:
                line = line.replace('COULCUTREPLACE', str(COUL_CUTOFF))
                
                
            f.write(line)
            f.write('\n')
            
            if 'neigh_modify' in line and write_groups:
                write_groups = False
                rigid_group = np.unique(groups)
                f.write('\n')
                f.write('group rig id ')
                for i in rigid_group:
                    f.write(str(int(i)) + ' ')
                f.write('\n')
                
                f.write('neigh_modify exclude molecule/intra rig')
                f.write('\n')

                
                f.write('group nonrigid subtract all rig')
                f.write('\n')
                
                f.write('\n#Nonbond Coeffs\n\n')
                for p in nb_params:
                    if len(p) == 6:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]} {p[4]} {p[5]}\n')
                    elif len(p) == 5:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]} {p[4]}\n')
                    elif len(p) == 4:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]}\n')
                    else:
                        raise ValueError('check nb params:', p)
                
                if fene:
                    f.write('\n#Bond Coeffs\n\n')
                    f.write(f'bond_coeff 1 fene 30.0 {1.5*graft_bead_size} 1.0 {1.0*graft_bead_size}\n')
                    if embed_distance == 0.63:
                        f.write('bond_coeff 2 harmonic 1000.0 0.37\n')
                    else:
                        f.write(f'bond_coeff 2 harmonic 30 {(graft_bead_size + endbead_size)/2}\n')

                    f.write('special_bonds fene\n')
                else:
                    f.write('\n#Bond Coeffs\n\n')
                    for p in b_params:
                        print()
                        f.write(f'bond_coeff {int(p[0])} {p[1]} {p[2]}\n')
    
    
    in_file = 'pgp_prod.in'
    with open(template_prod, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    
    #WRITE INPUT PROD SCRIPT
    write_groups = True
    nb_params = make_nb_params(graft_wca)
    
    with open(in_file,'w') as f:
        for line in lines:            
            if 'DIELREPLACE' in line:
                line = line.replace('DIELREPLACE', str(dielectric))
            if 'DEBYEREPLACE' in line:
                line = line.replace('DEBYEREPLACE', str(1/debye_length)) #take inverse of debye length here
            if 'COULCUTREPLACE' in line:
                line = line.replace('COULCUTREPLACE', str(COUL_CUTOFF))
            
            if line=='group centers type 1 4' and sticky_ends:
                line='group centers type 1 5'
            
            if fene and line == 'bond_style harmonic':
                line = 'bond_style hybrid harmonic fene'
                
            f.write(line)
            f.write('\n')
            
            if 'neigh_modify' in line and write_groups:
                write_groups = False
                
                f.write('\n')
                f.write('neigh_modify exclude molecule/intra rig')
                f.write('\n')

                
                f.write('\n#Nonbond Coeffs\n\n')
                for p in nb_params:
                    if len(p) == 6:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]} {p[4]} {p[5]}\n')
                    elif len(p) == 5:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]} {p[4]}')
                        f.write('\n')
                    elif len(p) == 4:
                        f.write(f'pair_coeff {int(p[0])} {int(p[1])} {p[2]} {p[3]}')
                        f.write('\n')
                    else:
                        raise ValueError('check nb params')
                        
                if fene:
                    f.write('\n#Bond Coeffs\n\n')
                    f.write(f'bond_coeff 1 fene 30.0 {1.5*graft_bead_size} 1.0 {1.0*graft_bead_size}\n')
                    if embed_distance == 0.63:
                        f.write('bond_coeff 2 harmonic 1000.0 0.37\n')
                    else:
                        f.write(f'bond_coeff 2 harmonic 30 {(graft_bead_size + endbead_size)/2}\n')
                    f.write('special_bonds fene\n')
                else:
                    f.write('\n#Bond Coeffs\n\n')
                    for p in b_params:
                        f.write(f'bond_coeff {int(p[0])} {p[1]} {p[2]}\n')
                    
                

if __name__ == '__main__':
    do_setup()
    
                    
        
        
