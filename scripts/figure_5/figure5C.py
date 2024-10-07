import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest
import os

contingency_dir = '../Dibyachintan_et_al_2024/datasets/contingent_mutations/'
input_dir = '../Dibyachintan_et_al_2024/datasets/contingent_fates/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_5/'
stat_dir ='../Dibyachintan_et_al_2024/statistical_tests/figure_5/figure_5C/'

rsa_dir = '../Dibyachintan_et_al_2024/rsa_dssp/'
df_rsa = pd.read_csv(rsa_dir + 'Myo3_rsa.dat', sep='\t')
invariant_sites = [35, 36, 38, 48, 49, 51, 54]

list_residue = []
for index, row in df_rsa.iterrows():
    rsa = row['rsa']
    if rsa <= 0.25:
        list_residue.append('B')
    else:
        list_residue.append('E')

df_rsa['type'] = list_residue
df_rsa = df_rsa[['resnum', 'aa', 'rsa', 'type']].copy()
df_rsa.columns = ['position', 'wt_aa', 'rsa', 'rsa_type']
df_rsa = df_rsa[(df_rsa['position'] > 0) & (df_rsa['position'] < 60)].copy()

df_all_mut = pd.read_excel(contingency_dir + 'Osh2.xlsx', engine='openpyxl')
df_all_mut = pd.merge(df_all_mut, df_rsa, on=['position'])
df_all_mut = df_all_mut[~df_all_mut['position'].isin(invariant_sites)].copy()

files = os.listdir(input_dir)

for file in files:
    df = pd.read_excel(input_dir + file)
    df = pd.merge(df, df_rsa, on=['position'])

    list_pval = []
    list_zscore = []

    list_type = ['B', 'E']
    count = []
    nobs = []

    df_t = df[df['fate'] == 'NF'].copy()
    df_c = df_all_mut.copy()

    for idx1 in range(len(list_type)):
        type1 = list_type[idx1]

        c1 = len(df_t[df_t['rsa_type'] == type1])
        n1 = len(df_c[df_c['rsa_type'] == type1])

        count.append(c1)
        nobs.append(n1)

    count = np.array(count)
    nobs = np.array(nobs)
    stat, pval = proportions_ztest(count, nobs)

    list_type.append('NF')
    list_pval.append(pval)
    list_zscore.append(stat)

    perc1 = (count[0]/nobs[0])
    perc2 = (count[1]/nobs[1])

    list_type = ['B', 'E']
    count = []
    nobs = []

    df_t = df[df['fate'] == 'SF'].copy()
    for idx1 in range(len(list_type)):
        type1 = list_type[idx1]

        c1 = len(df_t[df_t['rsa_type'] == type1])
        n1 = len(df_c[df_c['rsa_type'] == type1])

        count.append(c1)
        nobs.append(n1)

    count = np.array(count)
    nobs = np.array(nobs)
    stat, pval = proportions_ztest(count, nobs)

    list_pval.append(pval)
    list_zscore.append(stat)

    perc3 = (count[0]/nobs[0])
    perc4 = (count[1]/nobs[1])


    list_residue = ['buried', 'exposed', 'buried', 'exposed']
    list_fate = ['CNF', 'CNF', 'CSF', 'CSF']
    list_perc = [perc1, perc2, perc3, perc4]

    df_pval = pd.DataFrame()
    df_pval['fate'] = ['CNF', 'CSF']
    df_pval['pval'] = list_pval

    df_fate = pd.DataFrame()
    df_fate['residue'] = list_residue
    df_fate['fate'] = list_fate
    df_fate['score'] = list_perc


    fig, ax = plt.subplots(ncols=1, figsize=(7, 6))

    ax_color = 'black'

    palette = ["black", "snow"]

    sns.barplot(data=df_fate, x="fate", y="score", hue="residue",
                ax=ax, palette=palette, linewidth=2.5, edgecolor='black')

    ax.set_ylim([0, 0.15])
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['contingency\nin\nnonfunctionalization', 'contingency\nin\nsubfunctionalization'],
                       weight='bold')
    ax.set_yticks([0, 0.05, 0.10, 0.15])
    ax.set_yticklabels(['0.00', '0.05', '0.10', '0.15'], weight='bold')
    # ax.set_xticks(xticks)
    ax.tick_params(axis='x', direction='out', size=18, length=0, width=0, colors=ax_color, labelsize=18)
    ax.tick_params(axis='y', direction='out', size=18, length=8, width=4, colors=ax_color, labelsize=18)
    ax.set_ylabel('proportion of mutations', fontsize=24, color=ax_color, weight='bold')
    ax.set_xlabel('', fontsize=20, color=ax_color, weight='bold')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_edgecolor(ax_color)
    # ax.spines['bottom'].set_edgecolor(ax_color)
    ax.spines['bottom'].set_visible(False)
    ax.spines['bottom'].set_linewidth(4)
    ax.spines['left'].set_linewidth(4)
    ax.legend(loc='upper left', fontsize='x-large', ncol=1, labelspacing=0.15)

    fig.tight_layout()
    plt.savefig(figure_dir + file.replace('.xlsx', '.pdf'), dpi=600, bbox_inches='tight')
    plt.close()

    df_pval.to_excel(stat_dir + file, index=False, header=True,  engine='openpyxl')
