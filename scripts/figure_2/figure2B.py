import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_2/'

# plotting function
def plotting(x, y, figsize=(8, 8), figure=''):

    x_min = np.min(x)
    y_min = np.min(y)
    abs_min = min(x_min, y_min)

    x_max = np.max(x)
    y_max = np.max(y)
    abs_max = max(x_max, y_max)
    lims = [abs_min - 0.05, abs_max + 0.05]

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    ax.scatter(x, y, marker='o', color='#75648c', zorder=1, alpha=0.7)

    ax.plot(lims, lims, 'k--', alpha=0.9, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel('', fontsize=24, weight='bold')
    ax.set_ylabel('', fontsize=24, weight='bold')

    ax.tick_params(axis='both', which='major', labelsize=24, length=10, width=4, color='black')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_edgecolor('black')
    ax.spines['bottom'].set_edgecolor('black')
    ax.spines['left'].set_linewidth(4)
    ax.spines['bottom'].set_linewidth(4)
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(['0.0', '0.5', '1.0'])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_yticklabels(['0.0', '0.5', '1.0'])

    ax.set_xlabel('${{\Delta}F_{MYO3}^{Myo3}}$', fontsize=24)
    ax.set_ylabel('${{\Delta}F_{MYO5}^{Myo3}}$', fontsize=24)

    plt.savefig(figure_dir + figure + '.pdf', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':

    # files containing functional effect of genetic backgrounds being analyzed
    file1 = 'Osh2_MYO3_Myo3.dat'
    file2 = 'Osh2_MYO5_Myo3.dat'

    # read the functional effect file into dataframes
    df1 = pd.read_csv(input_dir + file1, sep='\t')
    df2 = pd.read_csv(input_dir + file2, sep='\t')

    # annotating mutations in both dataframes for comparison
    df1['mutated_reference'] = df1['position'].astype(str) + df1['mut_aa']
    df2['mutated_reference'] = df2['position'].astype(str) + df2['mut_aa']

    df1 = df1[['mutated_reference', 'dF']]
    df2 = df2[['mutated_reference', 'dF']]

    df = pd.merge(df1, df2, on='mutated_reference')

    x = df['dF_x']
    y = df['dF_y']

    plotting(x, y, figsize=(8, 8), figure='figure_2B')