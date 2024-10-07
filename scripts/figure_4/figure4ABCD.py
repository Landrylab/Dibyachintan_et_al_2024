import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_dir = '../Dibyachintan_et_al_2024/datasets/contingent_mutations/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_4/figure_4D/'

sorter = ['P', 'C', 'G', 'A', 'V', 'L', 'I', 'M', 'W', 'F', 'Y', 'H',
          'K', 'R', 'S', 'T', 'N', 'Q', 'D', 'E']

df1 = pd.read_excel(input_dir + 'Osh2' + '.xlsx', engine='openpyxl')

for i in range(1, 60):

    df = df1[df1['position'] == i].copy()
    df.mut_aa = df.mut_aa.astype("category")
    df.mut_aa = df.mut_aa.cat.set_categories(sorter)
    df = df.sort_values(by=["mut_aa"], ascending=False)

    ordered_df = df.copy()
    yticks = df['mut_aa'].tolist()
    my_range = range(1, len(df.index)+1)

    fontsize = 20

    fig, ax = plt.subplots(figsize=(6, 20))

    ax.hlines(y=my_range, xmin=ordered_df['dF_x'], xmax=ordered_df['dF_y'], color='black',
              alpha=1, zorder=0, lw=4)
    ax.scatter(ordered_df['dF_x'], my_range, color='#037ab1',
               alpha=1, label=r'$\{f}_{myo3}$', zorder=1, s=200)
    ax.scatter(ordered_df['dF_y'], my_range, color='#49997c',
               alpha=1, label=r'$\{f}_{myo5}$', zorder=1, s=200)
    ax.set_yticks(list(range(1, len(df.index)+1)))
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
    # ax.set_ylim(0.65, 3.25)
    plt.savefig(figure_dir + str(i) + '.pdf', bbox_inches='tight', dpi=600)
    plt.close()