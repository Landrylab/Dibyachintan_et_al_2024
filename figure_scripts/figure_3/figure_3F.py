import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import pandas as pd
from statistics import mean
from statistics import stdev
import matplotlib as mpl

mpl.rcParams.update({'errorbar.capsize': 6})
mpl.rc('font', family='Helvetica')

hfont = {'fontname': 'Helvetica'}

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/dms_functional_score/'
phylo_dir = parent_dir + 'datasets/phylogenetic_dataset/'

columns = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

def plot_distribution(df, threshold, width, color, figsize, title):

    plt.style.use('seaborn-white')
    # Data
    ppi = df['prey'].tolist()
    anc = np.array(df['ancestral'].tolist())
    homolog = np.array(df['homolog'].tolist())
    non_homolog = np.array(df['non_homolog'])

    anc_se = np.array(df['ancestral_se'].tolist())
    homolog_se = np.array(df['homolog_se'].tolist())
    non_homolog_se = np.array(df['non_homolog_se'].tolist())

    y = np.arange(3)

    fig, ax = plt.subplots(ncols=1, figsize=figsize)

    ax.barh(y + width, anc, width, xerr=[(0, 0, 0), anc_se], color='black', edgecolor="black")
    ax.barh(y, homolog, width, xerr=[(0, 0, 0), homolog_se], color='#939598', edgecolor="black")
    ax.barh(y - width, non_homolog, width, xerr=[(0, 0, 0), non_homolog_se], color='white', edgecolor="black")


    ax.set_xlim([0, 25])
    ax.set_yticks(y)
    ax.set_yticklabels(ppi, weight='bold')
    ax.tick_params(axis='x', direction='out', length=8, width=4, colors=color, labelsize=18)
    ax.tick_params(axis='y', direction='out', length=0, width=0, colors=color, labelsize=18)
    ax.set_xlabel('percentage of contingent mutations', fontsize=20, color=color, weight='bold')
#    ax.legend(["ancestral", "homologous", "non-homologous"], fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_edgecolor(color)
    ax.spines['bottom'].set_linewidth(4)

    fig.tight_layout()
    plt.show()
    return fig, ax



anc_df = pd.read_csv(phylo_dir + 'pre_duplication_ancestral_states.dat', sep='\t')
homolog_df = pd.read_csv(phylo_dir + 'non_ancestral_homologous_states.dat', sep='\t')
extant_df = pd.read_csv(phylo_dir + 'extant_paralogous_states.dat', sep='\t')

anc_list = anc_df['mutated_reference'].tolist()
homolog_list = homolog_df['mutated_reference'].tolist()
extant_list = extant_df['mutated_reference'].tolist()

anc_list = list(set(anc_list) - set(extant_list))

homolog_list = list(set(homolog_list) - set(anc_list))

homolog_list = list(set(homolog_list) - set(extant_list))

all_homolog_list = anc_list.copy()
all_homolog_list.extend(homolog_list)

all_homolog_list = list(set(all_homolog_list))


preys = ['Bbc1', 'Osh2', 'Pan1']

# for MYO3 locus baits are
baits = ['MYO3', 'MYO35']
title = 'MYO3 locus'

# for MYO5 locus baits are
# baits = ['MYO53', 'MYO5']
# title = 'MYO5 locus'

# contingent mutation if absolute value of ddF is greater than threshold
threshold = 0.2



if __name__ == '__main__':

    mut_ref_list = []

    counter = 0

    for prey in preys:
        df1 = pd.read_csv(input_dir + prey + '_' + baits[0] + '.dat', sep='\t')
        df2 = pd.read_csv(input_dir + prey + '_' + baits[1] + '.dat', sep='\t')
        df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
        df2 = df2[df2['sequence_type'] != 'synonymous'].copy()

        df1['mut_ref'] = df1['position'].astype(str) + df1['mut_aa']
        df2['mut_ref'] = df2['position'].astype(str) + df2['mut_aa']

        df = pd.merge(df1, df2, on='mut_ref')
        mut_ref = df['mut_ref'].tolist()

        if counter == 0:
            mut_ref_list.extend(mut_ref)
        else:
            mut_ref_list = list(set(mut_ref_list).intersection(set(mut_ref)))

        counter = counter + 1

    dict_plot = dict()

    dict_fisher_exact_anc = dict()
    dict_fisher_exact_hom = dict()
    dict_fisher_exact_anc_non = dict()

    anc_p_list = []
    homolog_p_list = []
    non_homolog_p_list = []

    se_anc_list = []
    se_homolog_list = []
    se_non_homolog_list = []
    for prey in preys:
        df1 = pd.read_csv(input_dir + prey + '_' + baits[0] + '.dat', sep='\t')
        df2 = pd.read_csv(input_dir + prey + '_' + baits[1] + '.dat', sep='\t')
        df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
        df2 = df2[df2['sequence_type'] != 'synonymous'].copy()

        df1['mut_ref'] = df1['position'].astype(str) + df1['mut_aa']
        df2['mut_ref'] = df2['position'].astype(str) + df2['mut_aa']

        df = pd.merge(df1, df2, on='mut_ref')
        mut_ref = df['mut_ref'].tolist()


        df['ddF'] = df['dF_y'] - df['dF_x']
        df['SED'] = np.sqrt(df['SE_dF_x'] ** 2 + df['SE_dF_y'] ** 2)
        df['1'] = df['dF1_y'] - df['dF1_x']
        df['2'] = df['dF2_y'] - df['dF1_x']
        df['3'] = df['dF3_y'] - df['dF1_x']
        df['4'] = df['dF1_y'] - df['dF2_x']
        df['5'] = df['dF2_y'] - df['dF2_x']
        df['6'] = df['dF3_y'] - df['dF2_x']
        df['7'] = df['dF1_y'] - df['dF3_x']
        df['8'] = df['dF2_y'] - df['dF3_x']
        df['9'] = df['dF3_y'] - df['dF3_x']


        anc_df = df[df['mut_ref'].isin(anc_list)].copy()
        homolog_df = df[df['mut_ref'].isin(homolog_list)].copy()
        non_homolog_df = df[~df['mut_ref'].isin(all_homolog_list)].copy()
        non_homolog_df = non_homolog_df[~non_homolog_df['mut_ref'].isin(extant_list)].copy()

        anc_df1 = anc_df[(anc_df['classification_x'] == 'binding') | (anc_df['classification_y'] == 'binding')].copy()
        homolog_df1 = homolog_df[(homolog_df['classification_x'] == 'binding') | (homolog_df['classification_y'] == 'binding')].copy()
        non_homolog_df1 = non_homolog_df[(non_homolog_df['classification_x'] == 'binding') | (non_homolog_df['classification_y'] == 'binding')].copy()

        anc_df1 = anc_df1[abs(anc_df1['ddF']) > 2 * anc_df1['SED']].copy()
        homolog_df1 = homolog_df1[abs(homolog_df1['ddF']) > 2 * homolog_df1['SED']].copy()
        non_homolog_df1 = non_homolog_df1[abs(non_homolog_df1['ddF']) > 2 * non_homolog_df1['SED']].copy()

        anc_p_list_x = []
        homolog_p_list_x = []
        non_homolog_p_list_x = []

        for col in columns:
            anc_df2 = anc_df1[abs(anc_df1[col]) > threshold].copy()
            homolog_df2 = homolog_df1[abs(homolog_df1[col]) > threshold].copy()
            non_homolog_df2 = non_homolog_df1[abs(non_homolog_df1[col]) > threshold].copy()

            perc_anc = round((len(anc_df2.index) / len(anc_df.index)) * 100, 2)
            perc_homolog = round((len(homolog_df2.index) / len(homolog_df.index)) * 100, 2)
            perc_non_homolog = round((len(non_homolog_df2.index) / len(non_homolog_df.index)) * 100, 2)

            anc_p_list_x.append(perc_anc)
            homolog_p_list_x.append(perc_homolog)
            non_homolog_p_list_x.append(perc_non_homolog)

        anc_p_list.append(mean(anc_p_list_x))
        homolog_p_list.append(mean(homolog_p_list_x))
        non_homolog_p_list.append(mean(non_homolog_p_list_x))

        anc_val = round((mean(anc_p_list_x) * len(anc_df.index)) / 100, 0)
        hom_val = round((mean(homolog_p_list_x) * len(homolog_df.index)) / 100, 0)
        non_val = round((mean(non_homolog_p_list_x) * len(non_homolog_df.index) / 100), 0)

        se_anc_list.append(stdev(anc_p_list_x))
        se_homolog_list.append(stdev(homolog_p_list_x))
        se_non_homolog_list.append(stdev(non_homolog_p_list_x))

        a = anc_val
        b = hom_val
        c = len(anc_df.index) - anc_val
        d = len(homolog_df.index) - hom_val

        res = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        dict_fisher_exact_anc[prey] = res[1]

        a = hom_val
        b = non_val
        c = len(homolog_df.index) - hom_val
        d = len(non_homolog_df.index) - non_val

        res = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        dict_fisher_exact_hom[prey] = res[1]

        a = anc_val
        b = non_val
        c = len(anc_df.index) - anc_val
        d = len(non_homolog_df.index) - non_val

        res = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        dict_fisher_exact_anc_non[prey] = res[1]

    df = pd.DataFrame()

    df['prey'] = preys
    df['ancestral'] = anc_p_list
    df['homolog'] = homolog_p_list
    df['non_homolog'] = non_homolog_p_list
    df['ancestral_se'] = se_anc_list
    df['homolog_se'] = se_homolog_list
    df['non_homolog_se'] = se_non_homolog_list


    plot_distribution(df, threshold, 0.2, '#3b444b', (6, 5), title)

    list_fe = [dict_fisher_exact_anc, dict_fisher_exact_anc_non, dict_fisher_exact_hom]
    df_fe = pd.DataFrame(list_fe)