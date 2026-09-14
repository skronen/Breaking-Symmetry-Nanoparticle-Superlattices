import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import plot_params


sav = True

savedir = 'paper_plots'

df_sa_5 = pd.read_csv('../results/diamseries_0.5.csv')
df_sa_8 = pd.read_csv('../results/diamseries_0.8.csv')
gdcol = r'Graft Density (Chains/nm$^2$)'
df_sa_5[gdcol] = 0.5
df_sa_8[gdcol] = 0.8

#%% plot for main series varying diameter for gd 0.8
df_asym = pd.read_csv('../results/diamseries_asym_0.8.csv')
f, ax = plt.subplots(figsize = (8, 5))
df_sym = df_sa_8.copy()
df_sym = df_sym[df_sym['diel']==15]
df_asym['meanfrac'] = (df_asym['frac1']+df_asym['frac2'])/2

df_asym_mean = df_asym.groupby('diam').mean()
df_asym_std = df_asym.groupby('diam').std()

df_sym_mean = df_sym.groupby('diam').mean()
df_sym_std = df_sym.groupby('diam').std()

means = np.concatenate([np.array(df_sym_mean['frac']), np.array(df_asym_mean['meanfrac'])]) *100
stds = np.concatenate([np.array(df_sym_std['frac']), np.array(df_asym_std['meanfrac'])]) *100

colors = np.array(['r', 'purple', 'b', 'g', 'g', 'goldenrod'])
col2 = ['limegreen', 'limegreen', 'gold']

frac1means = np.array(df_asym_mean['frac1'])*100
frac2means = np.array(df_asym_mean['frac2'])*100
frac1std = np.array(df_asym_std['frac1'])*100
frac2std = np.array(df_asym_std['frac2'])*100


labels = ['CsCl', 'NaCl 1', 'NaCl 2', 'NaCl 3', r'CaF$_2$', r'CaF$_2$', r'UH$_3$']
labels = ['NaCl', 'Th$_3$P$_4$', 'CsCl', r'CaF$_2$ (1)', r'CaF$_2$ (2)', r'UH$_3$']

for i in range(len(labels)):
    ax.errorbar(i, means[i], stds[i], color = colors[i], marker = 'o', 
            linestyle = 'none', capsize = 3, markersize = 7, capthick = 2, linewidth = 2)
    
    if i >= 3:
        ix = i-3
        ax.errorbar(i, frac1means[ix], frac1std[ix], color = col2[ix], marker = '^', 
                linestyle = 'none', capsize = 3, markersize = 7, capthick = 2, linewidth = 2)
        ax.errorbar(i, frac2means[ix], frac2std[ix], color = col2[ix], marker = '^', 
                linestyle = 'none', capsize = 3, markersize = 7, capthick = 2, linewidth = 2)
    
ax.set_ylabel('% Surface Area Forming Bonds')
ax.set_xticks(range(len(labels)), labels)
ax.set_xlim(-0.5, len(labels)-0.5)
ax.set_ylim(0, 35)
f.tight_layout()

if sav:
    f.savefig(os.path.join(savedir, 'asym_fracs.png'), dpi = 400)
    plt.close()
    
"""
compute asymmetric ratios
"""
val1 = 100/frac1means
val2 = 100/frac2means

vals = np.vstack((val1, val2))
maxs = np.max(vals, axis = 0)
mins = np.min(vals, axis = 0)

ratios = maxs/mins
avgs = np.mean(vals, axis = 0)

print(avgs, ratios)