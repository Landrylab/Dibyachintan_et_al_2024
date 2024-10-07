import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

parent_dir =  '../Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/functional_scores_chimera/'
figure_dir = parent_dir + 'figures/figure_6/'

list_colors_y = ['#ed1c24', '#231f20', '#231f20',
                 '#ed1c24', '#231f20', '#231f20', '#231f20',
                 '#ed1c24', '#ed1c24', '#231f20', '#ed1c24']

list_colors_x = ['#231f20', '#ed1c24', '#ed1c24', '#ed1c24']

def heatmap_contingency(data, title, yticklabels, xticklabels):

    fig, ax = plt.subplots(1, 1, figsize=(4, 9))

    fontsize = 18

    # vmin=-4 for Osh2, vmin=-3 for Bbc1
    g = sns.heatmap(data, cmap='bwr_r', linecolor='white', linewidth=5, center=0,
                    vmax=1, vmin=-3, square=False, cbar=True,
                    ax=ax,
                    cbar_kws={"orientation": "horizontal",
                              "shrink": 0.75, "pad": 0.08,
                              "aspect": 12})

    g.set_facecolor('grey')

    ax.set_yticklabels(yticklabels, weight='bold', fontsize=20)
    ax.set_xticklabels(xticklabels, weight='bold', fontsize=20)



    ax.tick_params(axis='y', which='major', rotation = 0,
                   labelsize=fontsize, length=0, width=4, color='#a7ba42')

    ax.tick_params(axis='x', which='major', rotation=0,
                   labeltop=True, labelbottom=False,
                   labelsize=fontsize, length=0, width=4, color='#a7ba42')



    ax.set_xlabel('', fontsize=fontsize, weight='bold')
    ax.set_ylabel('', fontsize=fontsize, weight='bold')

    cbar = ax.collections[0].colorbar
    for spine in cbar.ax.spines.values():  # show the colorbar frame again
        spine.set(visible=True, lw=.8, edgecolor='black')

    # Bbc1 colorscale
    cbar.set_ticks([-3, -2, -1, 0, 1])
    cbar.set_ticklabels(['-3', '-2', '-1', '0', '1'], fontsize=16, weight='bold')

    cbar.ax.tick_params(axis='x', direction='inout', colors='black', length=0)
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(2)

    labels = ax.get_yticklabels()
    ticks = ax.get_yticks()
    for label, tick in zip(labels, ticks):
        label.set_color(list_colors_y[int(tick-0.5)])
        label.set_fontname("Courier New")

    labels = ax.get_xticklabels()
    ticks = ax.get_xticks()
    for label, tick in zip(labels, ticks):
        label.set_color(list_colors_x[int(tick-0.5)])
        label.set_fontname("Courier New")

    plt.savefig(figure_dir + '6E.pdf', dpi=600, bbox_inches='tight')
    plt.show()

files = os.listdir(input_dir)

list_index = ['2M', '15P', '25V', '26I', '27Y',
              '29T', '31E', '39G', '56K', '58H', '59S']

yticklabels = ['k2M', 'S15p', 'I25v', 'v26I', 'F27y',
              'S29t', 'D31e', 'a39G', 't56K', 'Y58h', 'k59S']

xticklabels = ['wt', 'v26I', 'a39G', 't56K']

if __name__ == '__main__':

    for file in files:

        if 'MYO3_BBC1' not in file:
            continue

        df1 = pd.read_csv(input_dir + file, sep='\t')

        df1['mutated_reference_1'] = df1['position_1'].astype(str) + df1['AA_1']
        df1['mutated_reference_2'] = df1['position_2'].astype(str) + df1['AA_2']

        df1 = df1[['F', 'mutated_reference_1', 'mutated_reference_2', 'position_1', 'position_2']].copy()
        df2 = df1.copy()

        df2 = df2[['F', 'mutated_reference_2', 'mutated_reference_1', 'position_2', 'position_1']].copy()
        df2.columns = ['F', 'mutated_reference_1', 'mutated_reference_2', 'position_1', 'position_2']

        df_list = [df1, df2]

        df = pd.concat(df_list)

        df = df.sort_values(by=['position_1', 'position_2'], ascending=[True, True])

        df.loc[df.mutated_reference_1 == '0X', 'mutated_reference_1'] = 'WT'
        df.loc[df.mutated_reference_2 == '0X', 'mutated_reference_2'] = 'WT'

        list = [0, 26, 39, 56]

        df = df[~((df['mutated_reference_1'] == 'WT') & (df['mutated_reference_2'] == 'WT'))].copy()
        df = df[df['position_2'].isin(list)].copy()
        df = df[df['position_1'] != 0].copy()
        data = pd.pivot_table(df, values='F', index='mutated_reference_1', columns='mutated_reference_2')
        data = data[['WT', '26I', '39G', '56K']].copy()
        data = data.reindex(list_index)

        title = 'MYO3-Bbc1'

        heatmap_contingency(data, title, yticklabels, xticklabels)
