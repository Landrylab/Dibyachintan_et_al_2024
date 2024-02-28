import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/dms_variant_classification/'

list_exps = ['Aim21_MYO53', 'Aim21_MYO5',
             'Bbc1_MYO53', 'Bbc1_MYO5',
             'Pan1_MYO53', 'Pan1_MYO5',
             'Osh2_MYO53', 'Osh2_MYO5',
             'Arc18_MYO53', 'Arc18_MYO5']

list_exp1 = []
list_exp2 = []
list_r = []
list_rho = []
list_tau = []

for exp1 in list_exps:
    for exp2 in list_exps:

        df1 = pd.read_csv(input_dir + exp1 + '.dat', sep='\t')
        df2 = pd.read_csv(input_dir + exp2 + '.dat', sep='\t')


        df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
        df2 = df2[df2['sequence_type'] != 'synonymous'].copy()

        df1 = df1[['position', 'mut_aa', 'dF']].copy()
        df2 = df2[['position', 'mut_aa', 'dF']].copy()

        df1['reference'] = df1['position'].astype(str) + df1['mut_aa']
        df2['reference'] = df2['position'].astype(str) + df2['mut_aa']

        df3 = pd.merge(df1, df2, on='reference')
        df3 = df3.dropna()
        x = np.array(df3['dF_x'].tolist())
        y = np.array(df3['dF_y'].tolist())

        r, p1 = stats.pearsonr(x, y)
        rho, p2 = stats.spearmanr(x, y)
        tau, p3 = stats.kendalltau(x, y)

        exp1_list = exp1.split('_')
        exp2_list = exp2.split('_')

        if exp1_list[-1] == 'MYO5':
            list_exp1.append('MYO5    Myo5    ' + exp1_list[0])
        else:
            list_exp1.append('MYO5    Myo3    ' + exp1_list[0])

        if exp2_list[-1] == 'MYO53':
            list_exp2.append('MYO5    Myo3    ' + exp2_list[0])
        else:
            list_exp2.append('MYO5    Myo5    ' + exp2_list[0])




        list_r.append(r)
        list_rho.append(rho)
        list_tau.append(tau)


df = pd.DataFrame()
df['exp1'] = list_exp1
df['exp2'] = list_exp2
df['tau'] = list_tau
df['r'] = list_r
df['rho'] = list_rho

corr = 'tau'

# dms score correlation cluster heatmap code
df1 = df.copy()

sns.set(font_scale=1.4)
sns.set_style('white')

df_h = df1.copy()
df_heatmap_mean_ppi = pd.pivot_table(df_h, index='exp1', columns='exp2', values=corr)


colormap_sns = 'Purples'


res = sns.clustermap(df_heatmap_mean_ppi, mask=df_heatmap_mean_ppi.isna(),
                     vmin=0.5,
                     vmax=1.0,
                     cmap=colormap_sns,
                     linewidths=1, linecolor='WHITE',
                     cbar_kws={'label': ''}, figsize=(12, 12),
                     dendrogram_ratio=(.1, .10),
                     cbar_pos=(0.04, .9, .012, 0.08), col_cluster=True,
                     tree_kws={'linewidths': 5, 'zorder': 0},
                     metric='correlation', method='average')


res.ax_heatmap.set_xlabel('', fontsize=24, weight='bold', y=1.04)
res.ax_heatmap.set_ylabel('', fontsize=24, weight='bold', x=-0.04)

res.ax_heatmap.set_xticklabels([], weight='bold', fontsize=14,
                               rotation=90)
res.ax_heatmap.set_yticklabels(res.ax_heatmap.get_yticklabels(), weight='bold', fontsize=18)
res.ax_heatmap.tick_params(axis='y', direction='out', colors='black', length=0, width=0, labelsize=18)


colorbar = res.ax_heatmap.collections[0].colorbar
colorbar.set_ticks([0.5, 1.0])
colorbar.set_ticklabels(['0.5', '1.0'], weight='bold')

colorbar.ax.tick_params(axis='y', direction='out', colors='black', length=6, width=3, labelsize=0)
colorbar.ax.set_ylabel('', weight='bold', fontsize=18)
colorbar.ax.yaxis.set_label_position("left")
colorbar.outline.set_edgecolor('black')
colorbar.outline.set_linewidth(3)

plt.show()


