import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
import math

input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_3/figure_3B/'

myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
preys = ['Aim21', 'Arc18', 'Bbc1', 'Osh2', 'Pan1']
baits = ['MYO3_Myo3', 'MYO3_Myo5', 'MYO5_Myo3', 'MYO5_Myo5']


for bait in baits:

    counter = 0
    for prey in preys:
        df1 = pd.read_csv(input_dir + prey + '_' + bait + '.dat', sep='\t')

        df1['mutated_reference'] = df1['position'].astype(str) + df1['mut_aa']
    #    df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
    #    df1['dF'] = np.clip(df1['dF'].values, a_min=0.0, a_max=None)

        if counter == 0:
            df1['dF'] = np.clip(df1['dF'].values, a_min=0.0, a_max=1.0)
            df1 = df1[['mutated_reference', 'mut_aa', 'position', 'dF']].copy()
            df1.columns = ['mutated_reference', 'mut_aa', 'position', prey]
            df = df1.copy()
        else:
            df1['dF'] = np.clip(df1['dF'].values, a_min=0.0, a_max=1.0)
            df1 = df1[['mutated_reference', 'dF']].copy()
            df1.columns = ['mutated_reference', prey]
            df = pd.merge(df, df1, on='mutated_reference')

        counter = counter + 1


    groups = df.groupby('position')

    for position, df_t in groups:

        df_t = groups.get_group(position)
        df_h = df_t[['mut_aa', 'Aim21', 'Arc18', 'Bbc1', 'Osh2', 'Pan1']].copy()
        df_h.set_index('mut_aa', inplace=True)
        list_index = ['G', 'A', 'V', 'L', 'I', 'M', 'C', 'P', 'W', 'F', 'Y', 'H', 'K', 'R', 'S', 'T', 'N', 'Q', 'D', 'E']
        # list_index = ['G', 'A', 'V', 'L', 'I', 'M', 'C', 'P', 'W', 'F', 'Y', 'H', 'K', 'R', 'T', 'N', 'Q', 'D', 'E']
        df_heatmap_ppi = df_h.copy()
        df_heatmap_ppi = df_heatmap_ppi.reindex(list_index)

        fig, ax1 = plt.subplots(1, figsize=(4, 10))
        res = sns.heatmap(df_heatmap_ppi, mask=df_heatmap_ppi.isna(),
                          vmin=0.0,
                          vmax=1.0,
                          center=0.5, cmap='bwr_r',
                          linewidths=4, square=True, linecolor='white',
                          ax=ax1,
                          cbar=False)

        res.set_facecolor('grey')
        res.set_xlabel('position reference', fontsize=18, weight='bold')
        res.set_xlabel('', fontsize=18, weight='bold')
        res.set_ylabel('', fontsize=18, weight='bold')
        res.set_xticklabels(res.get_xmajorticklabels(), weight='bold', fontsize=14)
        res.set_yticklabels(res.get_ymajorticklabels(), weight='bold', fontsize=14)

        for tick in ax1.get_yticklabels():
            tick.set_fontname("Courier New")

        for tick in ax1.get_xticklabels():
            tick.set_fontname("Courier New")

        ax1.tick_params(axis='x', which='major', labelsize=18, length=0, width=0, color='black', rotation=90)
        ax1.tick_params(axis='y', which='major', labelsize=18, length=0, width=0, color='black', rotation=360)

        ax1.hlines([20], *ax1.get_xlim(), color='black')
        ax1.vlines([5], *ax1.get_ylim(), color='black')
        ax1.hlines([0], *ax1.get_xlim(), color='black')
        ax1.vlines([0], *ax1.get_ylim(), color='black')

        title = myo3_sh3_seq[position-1] + str(position)
        ax1.set_title(title, fontsize=20, weight='bold', y=1.01)
        plt.savefig(figure_dir + title + '_' + bait + '.pdf', dpi=600, bbox_inches='tight')
        plt.close()
