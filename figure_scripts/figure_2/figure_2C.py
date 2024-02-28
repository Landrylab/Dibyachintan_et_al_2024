import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/dms_variant_classification/'

myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
myo5_sh3_seq = 'PMFEAAYDFPGSGSPSELPLKKGDVIYITREEPSGWSLGKLLDGSKEGWVPTAYMKPHS'

df_median = pd.DataFrame()

preys = ['Osh2', 'Pan1', 'Arc18', 'Aim21', 'Bbc1']
# set sequence based on which SH3 domain is being plotted
seq = myo3_sh3_seq


for prey in preys:
    # Myo3 SH3 in MYO3 = _MYO3
    # Myo5 SH3 in MYO3 = _MYO35
    # Myo3 SH3 in MYO5 = _MYO53
    # Myo5 SH3 in MYO5 = _MYO5

    df1 = pd.read_csv(input_dir + prey + '_MYO3' + '.dat', sep='\t')

    df1 = df1[df1['sequence_type'] == 'non-synonymous'].copy()

    df3 = df1[['position', 'dF']].copy()

    df3['dF'] = np.clip(df3['dF'].values, a_min=0.0, a_max=1.0)

    df4 = pd.DataFrame(df3.groupby(['position'])['dF'].median())

    df_median[prey] = df4['dF']

    df_median['position'] = df_median.index

df_median = df_median[['position',
                       'Osh2', 'Pan1', 'Arc18', 'Aim21', 'Bbc1']].copy()

list_ref = []
for index, row in df_median.iterrows():
    pos = int(row['position'])
    if pos < 10:
        list_ref.append('00' + str(pos) + ' ' + seq[pos-1])
    else:
        list_ref.append('0' + str(pos) + ' ' + seq[pos - 1])


df_median['pos_ref'] = list_ref

df_median = df_median[['pos_ref',
                       'Osh2', 'Pan1', 'Arc18', 'Aim21', 'Bbc1']].copy()


df_median = df_median.set_index('pos_ref')

df_heatmap_ppi = df_median.transpose()
fig, ax = plt.subplots(1, 1, figsize=(30, 5))

colormap_sns = sns.diverging_palette(150, 275, s=90, l=40, as_cmap=True, center='light', sep=40)


res = sns.heatmap(df_heatmap_ppi, mask=df_heatmap_ppi.isna(),
                  vmin=0.0,
                  vmax=1.0,
                  center=0.5, cmap='bwr_r',
                  linewidths=5, square=False, linecolor='white',
                  ax=ax,
                  cbar=True)

colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0.0, 0.5, 1.0])
colorbar.set_ticklabels(['0.0', '0.5', '1.0'], fontsize=16, weight='bold')
colorbar.ax.tick_params(axis='y', direction='inout', colors='black', length=0)
colorbar.ax.set_ylabel(r'$\bf{\Delta}F_{median}$', weight='bold', fontsize=24)
colorbar.ax.yaxis.set_label_position("left")
colorbar.outline.set_edgecolor('black')
colorbar.outline.set_linewidth(1)

xlabels = df_median.index
res.set_facecolor('white')
ax.set_xlabel('', weight='bold', fontsize=20)
res.set_xticklabels(xlabels, weight='bold', rotation=90, fontsize=20, color='grey')
res.set_yticklabels(res.get_ymajorticklabels(), weight='heavy', rotation=0, fontsize=26)

chimera_list = [2, 15, 25, 26, 27, 29, 31, 39, 56, 58, 59]

for i in chimera_list:
    ax.get_xticklabels()[i-1].set_color("black")

ax.tick_params(axis='both', length=0)

ax.tick_params(bottom=True, top=False, left=True, right=False)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.show()