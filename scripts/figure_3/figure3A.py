import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt


input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_3/'

myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
myo5_sh3_seq = 'PMFEAAYDFPGSGSPSELPLKKGDVIYITREEPSGWSLGKLLDGSKEGWVPTAYMKPHS'

files = os.listdir(input_dir)

file = 'Bbc1_MYO3_Myo3.dat'

df1 = pd.read_csv(input_dir + file, sep='\t')

list_pos_ref = []
for index, row in df1.iterrows():
    if row['position'] < 10:
        list_pos_ref.append('00' + str(row['position']) + ' ' + row['wt_aa'])
    else:
        list_pos_ref.append('0' + str(row['position']) + ' ' + row['wt_aa'])

df1['position_reference'] = list_pos_ref

sns.set(font_scale=1.4)
sns.set_style('white')
df_h = df1[['position_reference', 'mut_aa', 'dF']].copy()

#
df_h['dF'] = np.clip(df_h['dF'].values, a_min=0.0, a_max=1)

# functional effects heatmap table
df_heatmap_ppi = pd.pivot_table(df_h, index='mut_aa', columns='position_reference', values='dF')
list_index = ['G', 'A', 'V', 'L', 'I', 'M', 'C', 'P', 'W', 'F', 'Y', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'D',
              'E']
df_heatmap_ppi = df_heatmap_ppi.reindex(list_index)

df_annotate = df1[['position_reference', 'mut_aa', 'sequence_type']].copy()
list_annot = []

# annotate synonymous mutations
for index, row in df_annotate.iterrows():
    if row['sequence_type'] == 'synonymous':
        list_annot.append('S')
    else:
        list_annot.append(' ')

df_annotate['annot'] = list_annot

# annotation heatmap table
df_annotate_ppi = pd.pivot_table(df_annotate, index='mut_aa', columns='position_reference', values='annot',
                                 aggfunc=np.sum)
df_annotate_ppi = df_annotate_ppi.reindex(list_index)


fig, ax1 = plt.subplots(1, figsize=(30, 10), dpi=300)
res = sns.heatmap(df_heatmap_ppi, mask=df_heatmap_ppi.isna(),
                  vmin=0.0,
                  vmax=1.0,
                  center=0.5, cmap='Purples_r',
                  linewidths=4, square=True, linecolor='white',
                  annot=df_annotate_ppi, fmt='', ax=ax1,
                  annot_kws={"fontsize": 14, "weight": 'bold'},
                  cbar=False)
#                      cbar_kws={'shrink': 0.2, 'aspect': 10, 'pad': 0.02})

# colorbar = ax1.collections[0].colorbar
# colorbar.set_ticks([0.0, 0.5, 1])
# colorbar.set_ticklabels(['0.0', '0.5', '1'])
# colorbar.ax.tick_params(axis='y', direction='inout', colors='black', length=8)
#
# colorbar.ax.set_ylabel('$\Delta$F', weight='bold')
# colorbar.ax.yaxis.set_label_position("left")
# colorbar.outline.set_edgecolor('black')
# colorbar.outline.set_linewidth(1)

res.set_facecolor('grey')
res.set_xlabel('position reference', fontsize=18, weight='bold')
res.set_xlabel('', fontsize=18, weight='bold')
res.set_ylabel('mutant AA', fontsize=18, weight='bold')
res.set_xticklabels(res.get_xmajorticklabels(), weight='bold', fontsize=14)
res.set_yticklabels(res.get_ymajorticklabels(), weight='bold', fontsize=14)

ax1.hlines([20], *ax1.get_xlim(), color='black')
ax1.vlines([59], *ax1.get_ylim(), color='black')
ax1.hlines([0], *ax1.get_xlim(), color='black')
ax1.vlines([0], *ax1.get_ylim(), color='black')

for tick in ax1.get_yticklabels():
    tick.set_fontname("Courier New")

for tick in ax1.get_xticklabels():
    tick.set_fontname("Courier New")

ax1.tick_params(axis='x', which='major', labelsize=18, length=0, width=0, color='black', rotation=90)
ax1.tick_params(axis='y', which='major', labelsize=18, length=0, width=0, color='black', rotation=360)

title = 'Bbc1 - MYO3 | Myo3 SH3'
ax1.set_title(title, fontsize=26, weight='bold', y=1.01)

plt.savefig(figure_dir + 'figure_3A.pdf', dpi=600, bbox_inches='tight')
plt.close()







