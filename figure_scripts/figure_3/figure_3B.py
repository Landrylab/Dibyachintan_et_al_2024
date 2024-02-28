import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.colors as mcolors
import matplotlib.colorbar as colorbar

import pandas as pd

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/statistical_dataset/functional_comparison_dataset/'


preys = ['Pan1', 'Osh2', 'Bbc1']

# change bait to MYO3 or MYO5
baits = ['MYO3']

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

def stacked_bar_chart(results, category_names):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    plt.style.use('seaborn-white')

    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)

    category_colors = [plt.cm.tab20c(0), plt.cm.tab20c(1), plt.cm.tab20c(2), plt.cm.tab20c(3),
                       plt.cm.Pastel1(8),
                       plt.cm.tab20c(11), plt.cm.tab20c(10), plt.cm.tab20c(9), plt.cm.tab20c(8)]

    cmap_col = category_colors
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.invert_yaxis()

    ax.xaxis.set_visible(True)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.9,
                        label=colname, color=color)

        r, g, b, _ = color
    #        text_color = 'white' if r * g * b < 0.5 else 'darkgrey'
    #        ax.bar_label(rects, label_type='center', color=text_color)
    # ax.legend(ncols=len(category_names), bbox_to_anchor=(0, 1),
    #           loc='lower left', fontsize='small')

    ax.spines[:].set_visible(False)

    fontsize = 36
    color = '#3b444b'
    ax.set_xlabel('percentage of functional mutations',
                   fontsize=fontsize, weight='bold', color=color)
    ax.set_xlabel('',
                  fontsize=fontsize, weight='bold', color=color)

    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(['0', '25', '50', '75', '100'])

    ax.tick_params(left=False, bottom=False)
    ax.tick_params(axis="x", direction="in", pad=-8, labelsize=fontsize, colors=color)
    ax.tick_params(axis="y", direction="out", pad=4, labelsize=fontsize, colors=color)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(preys, weight='bold', fontsize=fontsize)

    plt.show()


if __name__ == '__main__':

    mut_ref_list = []

    counter = 0

    for prey in preys:
        for bait in baits:

            df1 = pd.read_csv(input_dir + prey + ' - ' + bait + '.dat', sep='\t')
            mut_ref = df1['mutated_reference'].tolist()

            if counter == 0:
                mut_ref_list.extend(mut_ref)
            else:
                mut_ref_list = list(set(mut_ref_list).intersection(set(mut_ref)))

            counter = counter + 1


    dict_plot = dict()
    category_names = ['5', '10', '20', '25', 'ns', '25', '20', '10', '5']

    for prey in preys:
        for bait in baits:

            df1 = pd.read_csv(input_dir + prey + ' - ' + bait + '.dat', sep='\t')

            df1 = df1[df1['mutated_reference'].isin(mut_ref_list)].copy()


            pvalues = df1['pvalue'].tolist()
            adj_pvalues = p_adjust_bh(pvalues)
            df1['p_adjust_BH'] = adj_pvalues

            df3 = df1[df1['ddF'] < -0.1].copy()
            df5 = df1[df1['ddF'] > 0.1].copy()


            df_t = df5[df5['p_adjust_BH'] <= 0.05].copy()
            x5_1 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df5[df5['p_adjust_BH'] <= 0.10].copy()
            x5_2 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df5[df5['p_adjust_BH'] <= 0.20].copy()
            x5_3 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df5[df5['p_adjust_BH'] <= 0.25].copy()
            x5_4 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df3[df3['p_adjust_BH'] <= 0.05].copy()
            x3_1 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df3[df3['p_adjust_BH'] <= 0.10].copy()
            x3_2 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df3[df3['p_adjust_BH'] <= 0.20].copy()
            x3_3 = round((len(df_t.index) / len(df1.index)) * 100, 2)

            df_t = df3[df3['p_adjust_BH'] <= 0.25].copy()
            x3_4 = round((len(df_t.index) / len(df1.index)) * 100, 2)


            x_ns = 100.00 - x3_4 - x5_4
            x = np.array([x3_1, x3_2-x3_1, x3_3-x3_2, x3_4-x3_3,
                          x_ns,
                          x5_4-x5_3, x5_3-x5_2, x5_2-x5_1, x5_1])

            dict_plot[prey + '-' + bait] = x

    stacked_bar_chart(dict_plot, category_names)
