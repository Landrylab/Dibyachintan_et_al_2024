import os
import scipy.stats as st
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import math

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/chimera_functional_score/'

def plotting(vals, title, min, maximum):

    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    xs = []
    for i in range(2):
        xs.append(
            np.random.normal(i + 1, 0.04, vals[i].shape[0]))  # adds jitter to the data points - can be adjusted

    width = 4

    if 'MYO3' in title:
        ax_color = '#037ab1'
    else:
        ax_color = '#49997c'

    labels = ['', '']
    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels=labels, weight='bold')
    ax.set_yticks([min, 0.0])
    ax.set_xticklabels(labels=[str(min), '0.0'], weight='bold')

    ax.tick_params(axis='both', which='major', labelsize=14, length=6, width=width, color=ax_color)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_edgecolor(ax_color)
    ax.spines['bottom'].set_edgecolor(ax_color)
    ax.spines['left'].set_linewidth(width)
    ax.spines['bottom'].set_linewidth(width)


    ax.set_ylim([min-0.1, maximum+0.1])
    plt.axhline(y=0.0, color='#231f20', linestyle='--', linewidth=1)

    boxprops = dict(linestyle='-', linewidth=2, color=ax_color)
    flierprops = dict(marker='o', markersize=0,
                      linestyle='none')
    whiskerprops = dict(color=ax_color, linestyle='-', linewidth=2)
    capprops = dict(color=ax_color, linestyle='-', linewidth=2)
    medianprops = dict(linewidth=2, linestyle='-', color=ax_color)

    palette = ['#231f20', '#ed1c24']
    plt.boxplot(vals, labels=labels, notch=False, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops,
                flierprops=flierprops, medianprops=medianprops, showmeans=False, widths=[0.5, 0.5], zorder=0)

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=1.0, color=c, zorder=1, s=20)

    plt.show()


files = os.listdir(input_dir)
anc_states = ['2K', '15P', '25V', '26V', '27Y', '29T', '31E', '39A', '56T', '58H', '59K']

if __name__=='__main__':
    list_title = []
    list_pval = []
    list_stat = []
    for file in files:
        title = file.replace('_', '-').replace('.dat', '')
        df1 = pd.read_csv(input_dir + file, sep='\t')
        df1['mut_ref1'] = df1['position_1'].astype(str) + df1['AA_1']
        df1['mut_ref2'] = df1['position_2'].astype(str) + df1['AA_2']

        df_ancestral = df1[df1['mut_ref1'].isin(anc_states) & df1['mut_ref2'].isin(anc_states)].copy()
        df_divergent = df1[(~df1['mut_ref1'].isin(anc_states)) & (~df1['mut_ref2'].isin(anc_states))
                           & (~df1['mut_ref1'].isin(['0X'])) & (~df1['mut_ref2'].isin(['0X']))].copy()

        df_list = [df_ancestral, df_divergent]
        df = pd.concat(df_list)

        minimum = df['f'].min() / 0.2
        maximum = df['f'].max()
        min = round(math.floor(minimum) * (0.2), 1)

        anc_arr = np.array(df_ancestral['f'].tolist())
        div_arr = np.array(df_divergent['f'].tolist())
        vals = [anc_arr, div_arr]
        plotting(vals, title, min, maximum)

        stat, p_val = st.mannwhitneyu(anc_arr, div_arr, alternative='two-sided')

        list_title.append(title)
        list_stat.append(round(stat, 4))
        list_pval.append(round(p_val, 4))

    df = pd.DataFrame()
    df['exp'] = list_title
    df['stat'] = list_stat
    df['pval'] = list_pval

