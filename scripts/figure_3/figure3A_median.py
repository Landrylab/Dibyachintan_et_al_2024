import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_3/'

myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
myo5_sh3_seq = 'PMFEAAYDFPGSGSPSELPLKKGDVIYITREEPSGWSLGKLLDGSKEGWVPTAYMKPHS'

df_median = pd.DataFrame()

preys = ['Bbc1']
# set sequence based on which SH3 domain is being plotted
seq = myo3_sh3_seq


for prey in preys:

    df1 = pd.read_csv(input_dir + prey + '_MYO3_Myo3' + '.dat', sep='\t')

    df1 = df1[df1['sequence_type'] == 'non-synonymous'].copy()

    df3 = df1[['position', 'dF']].copy()

    df3['dF'] = np.clip(df3['dF'].values, a_min=0.0, a_max=1.0)

    df_m = pd.DataFrame(df3.groupby(['position'])['dF'].median())

    df_m['position'] = df_m.index

df_median = df_m[['position', 'dF']].copy()
df_median.columns = ['position', 'Bbc1']

list_ref = []
for index, row in df_median.iterrows():
    pos = int(row['position'])
    if pos < 10:
        list_ref.append('00' + str(pos) + ' ' + seq[pos-1])
    else:
        list_ref.append('0' + str(pos) + ' ' + seq[pos - 1])


df_median['pos_ref'] = list_ref

df_median = df_median[['position', 'Bbc1']].copy()
df_median.columns = ['position', 'Bbc1']

df_median = df_median.set_index('position')

df_heatmap_ppi = df_median.transpose()
fig, ax = plt.subplots(1, 1, figsize=(30, 1.4))



res = sns.heatmap(df_heatmap_ppi, mask=df_heatmap_ppi.isna(),
                  vmin=0.0,
                  vmax=1.0,
                  center=0.5, cmap='Purples_r',
                  linewidths=5, square=True, linecolor='white',
                  ax=ax,
                  cbar=False)

# colorbar = ax.collections[0].colorbar
# colorbar.set_ticks([0.0, 0.5, 1.0])
# colorbar.set_ticklabels(['0.0', '0.5', '1.0'], fontsize=16, weight='bold')
# colorbar.ax.tick_params(axis='y', direction='inout', colors='black', length=0)
# colorbar.ax.set_ylabel(r'$\bf{\Delta}F_{median}$', weight='bold', fontsize=24)
# colorbar.ax.yaxis.set_label_position("left")
# colorbar.outline.set_edgecolor('black')
# colorbar.outline.set_linewidth(1)

rsa_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_etal_2024_NatComms/rsa_dssp/'
df_rsa = pd.read_csv(rsa_dir + 'Myo3_rsa.dat', sep='\t')

list_residue = []
for index, row in df_rsa.iterrows():
    rsa = row['rsa']
    if rsa <= 0.04:
        list_residue.append('B')
    elif rsa <= 0.25:
        list_residue.append('b')
    elif rsa <= 0.50:
        list_residue.append('e')
    else:
        list_residue.append('E')

df_rsa['type'] = list_residue
df_rsa = df_rsa[['resnum', 'aa', 'rsa', 'type']].copy()
df_rsa.columns = ['position', 'wt_aa', 'rsa', 'type']

df_rsa = df_rsa[df_rsa['position'] >= 1].copy()
df_rsa = df_rsa[df_rsa['position'] <= 59].copy()
list_1 = df_rsa['type'].tolist()

rsa_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_etal_2024_NatComms/evo_rate4site/'
df_rsa = pd.read_csv(rsa_dir + 'Myo3.dat', delim_whitespace=True)

list_residue = []
for index, row in df_rsa.iterrows():
    evo = float(row['SCORE'])
    if evo <= -1.0:
        list_residue.append('I')
    elif evo <= 0.0:
        list_residue.append('S')
    else:
        list_residue.append('F')

df_rsa['type'] = list_residue
df_rsa = df_rsa[['POS', 'SEQ', 'SCORE', 'type']].copy()
df_rsa.columns = ['position', 'wt_aa', 'evo_rate', 'type']
list_2 = df_rsa['type'].tolist()

xticklabels = []
for idx in range(len(list_1)):
    xticklabels.append(list_1[idx] + '\n' + list_2[idx])


res.set_facecolor('white')
ax.set_xlabel('', weight='bold', fontsize=20)
res.set_xticklabels(xticklabels, weight='bold', rotation=0, fontsize=20, color='grey')
res.set_yticklabels([], weight='heavy', rotation=0, fontsize=26)
for tick in ax.get_xticklabels():
    tick.set_fontname("Courier New")
# chimera_list = [2, 15, 25, 26, 27, 29, 31, 39, 56, 58, 59]
#
# for i in chimera_list:
#     ax.get_xticklabels()[i-1].set_color("black")

ax.tick_params(axis='both', length=0)

ax.tick_params(bottom=True, top=False, left=True, right=False)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
plt.savefig(figure_dir + 'median_Bbc1' + '.pdf', dpi=600, bbox_inches='tight')
plt.close()
