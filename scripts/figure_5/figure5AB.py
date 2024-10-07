import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_5/figure_5A/'
preys = ['Pan1', 'Osh2', 'Bbc1']

mutants = ['A', 'C', 'D', 'E', 'F', 'G',
           'H', 'I', 'K', 'L', 'M', 'N',
           'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']




for i in range(1, 60):
    for mut in mutants:
        mutant = str(i) + mut

        df3_list = []
        df5_list = []

        for prey in preys:
            df33 = pd.read_csv(input_dir + prey + '_MYO3_Myo3.dat', sep='\t')
            df33['mutated_reference'] = df33['position'].astype(str) + df33['mut_aa']
            df33 = df33[df33['mutated_reference'] == mutant].copy()
            df33['myosin'] = ['MYO3']
            df33['prey'] = [prey]

            df35 = pd.read_csv(input_dir + prey + '_MYO3_Myo5.dat', sep='\t')
            df35['mutated_reference'] = df35['position'].astype(str) + df35['mut_aa']
            df35 = df35[df35['mutated_reference'] == mutant].copy()
            df35['myosin'] = ['MYO35']
            df35['prey'] = [prey]


            df3_list.append(df33)
            df5_list.append(df35)


        df33 = pd.concat(df3_list)

        df35 = pd.concat(df5_list)


        df3 = pd.merge(df33, df35, on=['mutated_reference', 'prey'])



        df3 = df3[['mutated_reference', 'prey', 'dF_x', 'dF_y']]


        df = df3


        ordered_df = df.sort_values(by='prey')
        my_range = range(1, len(df.index)+1)
        yticks = ordered_df['prey'].tolist()

        if len(yticks) < 3:
            continue
        fontsize = 28
        # The horizontal plot is made using the hline function
        fig, ax = plt.subplots(figsize=(6, 6))
        #d09c2e myo3
        #ae3b25 myo5
        # ax.hlines(y=my_range, xmin=ordered_df['mean_f_syn_stop_x'], xmax=ordered_df['mean_f_syn_stop_y'], color='#d09c2e',
        #           alpha=1, zorder=0, lw=8)
        ax.hlines(y=my_range, xmin=ordered_df['dF_x'], xmax=ordered_df['dF_y'], color='black',
                  alpha=1, zorder=0, lw=8)
        ax.scatter(ordered_df['dF_x'], my_range, color='#037ab1',
                   alpha=1, label=r'$\{f}_{myo3}$', zorder=1, s=400)
        ax.scatter(ordered_df['dF_y'], my_range, color='#49997c',
                   alpha=1, label=r'$\{f}_{myo5}$', zorder=1, s=400)
        ax.set_yticks([1, 2, 3])
        ax.set_yticklabels(labels=yticks, weight='bold')
        ax.set_xticks([0, 0.4, 0.8, 1.2])
        ax.set_xticklabels([0, 0.4, 0.8, 1.2], weight='bold')
        ax.tick_params(axis='both', which='major', labelsize=fontsize, length=10, width=4, color='black')
        ax.set_xlabel('$\Delta{F}$', fontsize=fontsize, weight='bold')
        #ax.set_xlabel('', fontsize=fontsize, weight='bold')
        #ax.set_xlabel('', fontsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_edgecolor('black')
        ax.spines['bottom'].set_edgecolor('black')
        ax.spines['left'].set_linewidth(4)
        ax.spines['bottom'].set_linewidth(4)
        ax.set_ylim(0.65, 3.25)
        # ax.axvline(1, ls='--', color='black', alpha=0.5, lw=4)
        # ax.axvline(0.75, ls='--', color='black', alpha=0.5, lw=4)
        ax.set_title(mutant, fontsize=fontsize, weight='bold', y=1.05)
        #plt.savefig(illustrator_dir + 'figure_3.pdf', bbox_inches='tight', dpi=300)
        plt.savefig(figure_dir + mutant + '.pdf', bbox_inches='tight', dpi=600)
        plt.close()
        # plt.show()

