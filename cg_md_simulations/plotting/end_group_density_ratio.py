    
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
                  scale = 0.75, linestyles = '', errwidth = 2, capsize = 0.1)#marker = 'o',)
plt.xlabel('Nanoparticle Diameter, nm')
plt.ylabel('Density Ratio (Contact/Isolated)')
plt.tight_layout()

if sav:
    plt.savefig(os.path.join(savedir, 'md_endgroup_densities_gd0.8.png'), dpi = 400)
    plt.close()
    