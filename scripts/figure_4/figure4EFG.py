import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as stats

input_dir = '../Dibyachintan_et_al_2024/datasets/contingent_mutations/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_4/figure_4A/'

invariant_sites = [35, 36, 38, 48, 49, 51, 54]
prey = 'Pan1'
def Rhist(x, y, bins=None, xlab='', prey = '', color='w', edgecolor='k', figsize=(8, 6), offset=5):
    """
    Makes histograms that look like R
    Inputs:
    - x: a numpy array or pandas series
    - bins: number of bins, default (None) is mpl default
    - xlab: text label for x axis, default '' (empty)
    - savename: full name and path of saved figure,
      if '' (default) nothing saved
    - color: fill color of bars, default 'w' (white)
    - edgecolor: outline color of bars, default 'k' (black)
    - figsize: width, heighth of figure in inches (default 8x6)
    - offset: how far to separate axis, default=5
    """
    plt.style.use('seaborn-ticks')

    def adjust_spines(ax, spines, offset):
        """
        This is mostly from
        https://matplotlib.org/examples/pylab_examples/spine_placement_demo.html
        """
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward', offset))  # outward by offset points
                # spine.set_bounds(0, 1)
            else:
                spine.set_color('none')  # don't draw spine

        # turn off ticks where there is no spine
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            # no yaxis ticks
            ax.yaxis.set_ticks([])

        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            # no xaxis ticks
            ax.xaxis.set_ticks([])

    ax_color = 'black'
    ax_width = 4
    tick_width = 4
    tick_length = 4
    tick_size = 8
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(x, bins=bins, color='gray', edgecolor=edgecolor)
    ax.hist(y, bins=bins, color=color, edgecolor=edgecolor)
    adjust_spines(ax, ['left', 'bottom'], offset)
    ax.set_xlabel(xlab, fontsize=20)
    ax.set_ylabel('number of mutations', fontsize=20)

    ax.tick_params(axis='x', direction='out', size=tick_size, length=tick_length, width=tick_width, colors=ax_color, labelsize=18)
    ax.tick_params(axis='y', direction='out', size=tick_size, length=tick_length, width=tick_width, colors=ax_color, labelsize=18)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_edgecolor(ax_color)
    ax.spines['bottom'].set_edgecolor(ax_color)

    ax.spines['bottom'].set_linewidth(ax_width)
    ax.spines['left'].set_linewidth(ax_width)

    ax.set_xlim([-0.6, 0.6])
    ax.set_ylim([0, 300])

    plt.tight_layout()
    plt.savefig(figure_dir + prey + '.pdf', dpi=600, bbox_inches='tight')
    plt.show()



if __name__ == '__main__':
    df1 = pd.read_excel(input_dir + prey + '.xlsx', engine='openpyxl')
    df1 = df1[~df1['position'].isin(invariant_sites)].copy()

    df_x = df1.copy()
    df_y = df1[df1['contingent'] == 'N'].copy()

    x = np.array(df_x['ddF'].tolist())
    y = np.array(df_y['ddF'].tolist())

    bins = np.linspace(-0.6, 0.6, num=21)

    Rhist(x, y, bins=bins, xlab='${\Delta\Delta}F$', prey=prey)

