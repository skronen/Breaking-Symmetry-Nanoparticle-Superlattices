    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import plot_params


sav = True
savedir = 'paper_plots'

df_density = pd.read_csv('../results/density_ratios.csv')
    
df = df_density
df_plot = df[df['diel']==15.0]
df_plot = df_plot[df_plot['sig']==0.8]
hue = 'Graft Density\n' + r'(chains/nm$^2$)'
df_plot[hue] = df_plot['sig']


plt.figure()
sns.pointplot(df_plot, x = 'd', y = 'ratio', palette = ['k'], errorbar = 'sd',
                  scale = 0.75, linestyles = '', errwidth = 2, capsize = 0.1)
plt.xlabel('Nanoparticle Diameter, nm')
plt.ylabel('Density Ratio (Contact/Isolated)')

map = {11.1: 0, 16.8: 1, 25.5: 2}
for d, ratio in zip(df_plot['d'], df_plot['ratio']):
    plt.scatter(map[d]+0.1, ratio, color = 'gray', alpha = 0.7, s = 15)

plt.tight_layout()

if sav:
    plt.savefig(os.path.join(savedir, 'md_endgroup_densities_gd0.8.png'), dpi = 400)
    plt.close()
    