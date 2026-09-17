import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import plot_params


sav = True

savedir = 'paper_plots'

df_sa_5 = pd.read_csv('../results/diamseries_gd0.5.csv')
df_sa_8 = pd.read_csv('../results/diamseries_gd0.8.csv')
gdcol = r'Graft Density (Chains/nm$^2$)'
df_sa_5[gdcol] = 0.5
df_sa_8[gdcol] = 0.8

#%% plot for main series varying diameter for gd 0.8

plt.figure()
df = df_sa_8
df_plot = df[df['diel']==15]
df_plot['diam'] = df_plot['diam']/10
df_plot['frac'] = df_plot['frac']*100
sns.barplot(df_plot, x = 'diam', y = 'frac', palette = ['r', 'purple', 'b'],
            errorbar = 'sd', capsize=.1, errwidth = 1.5, errcolor = 'k', width = 0.6)
#plt.legend()
plt.xlabel('Nanoparticle Diameter, nm')
plt.ylabel('% Surface Area Forming Bonds')
plt.ylim([0,24])

map = {11.1: 0, 16.8: 1, 25.5: 2}
for d, ratio in zip(df_plot['diam'], df_plot['frac']):
    plt.scatter(map[d]+0.1, ratio, color = 'gray', alpha = 0.7, s = 15)

plt.tight_layout()


if sav:
    plt.savefig(os.path.join(savedir, 'md_valencies_diamseries.png'), dpi = 400)
    plt.close()
    
#%% compare gd 0.5 to gd 0.8

plt.figure()
df = df_sa_8
df = pd.concat((df, df_sa_5), axis = 0)
df_plot = df[df['diel']==15]
df_plot['diam'] = df_plot['diam']/10
df_plot['frac'] = df_plot['frac']*100

sns.barplot(df_plot, x = 'diam', y = 'frac', hue = gdcol, palette = ['maroon', 'red'],
            errorbar = 'sd', capsize=.1, errwidth = 1.5, errcolor = 'k', width = 0.6)


for i,hue in enumerate(np.unique(df_plot[gdcol])):
    if i == 0: shift = -0.05
    else: shift = 0.25
    
    df_hue = df_plot[df_plot[gdcol] == hue]
    map = {11.1: 0, 16.8: 1, 25.5: 2}
    for d, ratio in zip(df_hue['diam'], df_hue['frac']):
        plt.scatter(map[d]+shift, ratio, color = 'gray', alpha = 0.7, s = 15)

plt.xlabel('Nanoparticle Diameter, nm')
plt.ylabel('% Surface Area Forming Bonds')


plt.tight_layout()
if sav:
    plt.savefig(os.path.join(savedir, 'md_valencies_varygd.png'), dpi = 400)
    plt.close()
    
#%% compare results across dielectric constants

f, ax = plt.subplots()
df_plot = df_sa_8.copy()

hue = 'Nanoparticle Diameter, nm'
df_plot[hue] = df_plot['diam']/10
df_plot['frac'] = df_plot['frac']*100


for ms, diam, col in zip(['^', 'o', 's'], np.unique(df_plot['diam']), ['r','purple', 'blue']):
    df= df_plot.copy()
    df = df[df['diam']==diam]
    with plt.rc_context({'lines.linewidth': 0.7}):
        sns.pointplot(df, x = 'diel', y = 'frac', hue = hue, palette = [col], errorbar = 'sd',
                  scale = 1, errwidth = 2, capsize = 0.1, ax = ax, markers = [ms])

for coll in ax.collections:
    if isinstance(coll, matplotlib.collections.PathCollection):
        coll.set_sizes([50])   

for line in ax.lines:
    line.set_linewidth(1.75) 

colors = ['firebrick', 'indigo', 'navy']
for i,hue in enumerate(np.unique(df_plot['diam'])):
    shift = 0.1
    
    df_hue = df_plot[df_plot['diam'] == hue]
    map = {10: 0, 15: 1, 22.5: 2}
    for d, ratio in zip(df_hue['diel'], df_hue['frac']):
        plt.scatter(map[d]+shift, ratio, c = colors[i], alpha = 0.7, s = 15)

ax.set_ylabel('% Surface Area Forming Bonds')
ax.set_xlabel('Dielectric Constant')

leg = ax.get_legend()
if leg is not None:
    for handle in leg.legendHandles:

        # Scatter-style legend marker
        if isinstance(handle, matplotlib.collections.PathCollection):
            handle.set_sizes([50])   # area in points^2

        # Line-style legend marker
        elif isinstance(handle, matplotlib.lines.Line2D):
            handle.set_markersize(np.sqrt(50))
            handle.set_linewidth(3)
            
f.tight_layout()
if sav:
    f.savefig(os.path.join(savedir, 'percs_vary_diel.png'), dpi = 400)
    plt.close()
    