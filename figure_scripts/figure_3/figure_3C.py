import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
variant_dir = parent_dir + 'datasets/statistical_dataset/functional_comparison_dataset/'
null_dir = parent_dir + 'datasets/statistical_dataset/non_functional_dataset/'

all_sites = []

# edit this for getting the three ppi plots
ppi ='Bbc1'
# ppi ='Pan1'
# ppi ='Osh2'

list_prey = ['Bbc1', 'Pan1', 'Osh2']

for prey in list_prey:
    list_sites = []

    file = prey + ' - MYO3.dat'

    df1 = pd.read_csv(variant_dir + file, sep='\t')
    df2 = pd.read_csv(null_dir + file, sep='\t')

    df1 = df1[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
               'ddF', 'SE_dF_3', 'SE_dF_5']].copy()

    df1['SED'] = np.sqrt(df1['SE_dF_3']**2 + df1['SE_dF_5']**2)

    df1 = df1[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
               'ddF', 'SED']].copy()
    df1.columns = ['position', 'mutated_aa', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
                   'ddF', 'SED']

    df1['type'] = ['FDR']*len(df1.index)

    df2 = df2[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5']].copy()

    df2['ddF'] = [0.0]*len(df2.index)
    df2['SED'] = [0.0]*len(df2.index)

    df2.columns = ['position', 'mutated_aa', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
                   'ddF', 'SED']

    df2['type'] = ['NF']*len(df2.index)

    df_list = [df1, df2]

    df = pd.concat(df_list)

    df_syn = df[(df['sequence_type_3'] == 'synonymous')
                & (df['sequence_type_5'] == 'synonymous')].copy()

    df_temp = df.copy()

    #remove amino acids coded by a single codon
    df_temp = df_temp[df_temp['mutated_aa'] != 'M'].copy()
    df_temp = df_temp[df_temp['mutated_aa'] != 'W'].copy()

    list_syn = df_syn['mutated_reference'].tolist()

    df_check = df_temp[abs(df_temp['ddF']) >= 2*df_temp['SED']].copy()

    df_check = df_check[abs(df_check['ddF']) >= 0.2].copy()
    df_check = df_check[~df_check['mutated_reference'].isin(list_syn)].copy()

    pos_groups = df_check.groupby('position')

    for pos, df_p in pos_groups:
        if len(df_p.index) > 2:
            list_sites.append(pos)
        else:
            pass

    set_sites = set(list_sites)

    all_sites = list(set_sites.union(set(all_sites)))

list_sites = []


file_name = ppi + ' - MYO3.dat'

df1 = pd.read_csv(variant_dir + file_name, sep='\t')
df2 = pd.read_csv(null_dir + file_name, sep='\t')

df1 = df1[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
           'ddF', 'SE_dF_3', 'SE_dF_5']].copy()

df1['SED'] = np.sqrt(df1['SE_dF_3']**2 + df1['SE_dF_5']**2)

df1 = df1[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
           'ddF', 'SED']].copy()
df1.columns = ['position', 'mutated_aa', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
               'ddF', 'SED']

df1['type'] = ['FDR']*len(df1.index)

df2 = df2[['position_3', 'mut_aa_3', 'mutated_reference', 'sequence_type_3', 'sequence_type_5']].copy()

df2['ddF'] = [0.0]*len(df2.index)
df2['SED'] = [0.0]*len(df2.index)

df2.columns = ['position', 'mutated_aa', 'mutated_reference', 'sequence_type_3', 'sequence_type_5',
               'ddF', 'SED']

df2['type'] = ['NF']*len(df2.index)

df_list = [df1, df2]

df = pd.concat(df_list)


df = df[df['position'].isin(all_sites)]

df['ddF'] = np.clip(df['ddF'].values, a_min=-0.5, a_max=0.5)

jitter = 0.05



pos_groups = df.groupby('position')

fig, ax = plt.subplots(figsize=(5, 10))

counter = 0
for pos, df_plot in pos_groups:
    df_x_jitter = np.random.normal(loc=0, scale=jitter, size=len(df_plot.index)) + counter
    colors = []
    alphas = []
    for index, row in df_plot.iterrows():
        if (row['ddF'] > 0.2) and (abs(row['ddF']) > 2*row['SED']):
            colors.append('#49997c')
            alphas.append(0.9)
        elif (row['ddF'] < -0.2) and (abs(row['ddF']) > 2*row['SED']):
            colors.append('#0097c3')
            alphas.append(0.9)
        else:
            colors.append('grey')
            alphas.append(0.25)

    values = df_plot['ddF'].tolist()
    ax.scatter(values, df_x_jitter, s=60, marker='o', alpha=alphas, zorder=1, c=colors)
    counter += 1


ax.set_yticks(range(len(all_sites)))
ax.set_yticklabels(all_sites)
ax.set_ylim(-0.5, len(all_sites) - 0.5)

ax_color = '#d09c2e'
ax_color = 'black'
fontsize = 20
ax.set_xlim(-0.55, 0.55)

ax.set_xticks([-0.50, -0.25, 0.00, 0.25, 0.5])
ax.set_xticklabels(['-0.50', '-0.25', '0.00', '0.25', '0.50'], weight='bold')
# ax.set_xticks(xticks)
ax.tick_params(axis='x', direction='out', length=8, width=4, colors=ax_color, labelsize=18)
ax.tick_params(axis='y', direction='out', length=0, width=0, colors=ax_color, labelsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_edgecolor(ax_color)
ax.spines['bottom'].set_edgecolor(ax_color)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.tick_params(axis='both', which='major', labelsize=fontsize, length=10, width=4, color=ax_color)
ax.axvline(0.0, ls='--', color=ax_color, alpha=1, lw=2)

ax.spines['left'].set_bounds(0, len(all_sites)-1)
ax.spines['bottom'].set_bounds(-0.5, 0.5)

plt.show()
