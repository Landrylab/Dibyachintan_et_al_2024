import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
phylo_dir = parent_dir + 'datasets/phylogenetic_dataset/'

posterior_prob_dir = phylo_dir + 'ancestral_node_posterior_prob/'
ancestral_df = pd.read_csv(phylo_dir + 'pre_duplication_ancestral_states.dat', sep='\t')

files = os.listdir(posterior_prob_dir)
df_list = []

for file in files:
    df1 = pd.read_csv(posterior_prob_dir + file, sep='\t')
    df1['mutated_reference'] = df1['paralog_position'].astype(str) + df1['ancestral_state']
    df2 = df1[['mutated_reference', 'prob']].copy()
    df_list.append(df2)

df = pd.concat(df_list)

df = df.sort_values('prob', ascending=False).drop_duplicates('mutated_reference').sort_index()


df_final = pd.merge(df, ancestral_df, on='mutated_reference')

x = np.array(df_final['prob'].tolist())

fig, ax = plt.subplots(figsize=(5, 4))

ax_color = 'black'
fontsize = 20

ax.hist(x, bins=20, color='lightgrey', edgecolor='black')
ax.set_xlabel('maximum posterior probability \n along evolutionary trajectory', fontsize=fontsize, weight='bold')
ax.set_ylabel('ancestral states', fontsize=fontsize, weight='bold')



ax.tick_params(axis='x', direction='out', length=8, width=4, colors=ax_color, labelsize=18)
ax.tick_params(axis='y', direction='out', length=0, width=0, colors=ax_color, labelsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_edgecolor(ax_color)
ax.spines['left'].set_edgecolor(ax_color)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.tick_params(axis='both', which='major', labelsize=fontsize, length=10, width=4, color=ax_color)
ax.axvline(np.median(x), ls='--', color=ax_color, alpha=1, lw=2)
ax.set_ylim([-1, 61])
ax.set_xlim([-0.06, 1.06])
ax.spines['left'].set_bounds(0, 60)
ax.spines['bottom'].set_bounds(0, 1.0)
plt.tight_layout()
plt.show()