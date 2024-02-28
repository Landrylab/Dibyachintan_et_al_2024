import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import pandas as pd
import matplotlib as mpl

parent_dir = '/Users/sohamdibyachintan/Desktop/projects/Dibyachintan_et_al_2024/'
variant_dir = parent_dir + 'datasets/statistical_dataset/comparison_dataset/'
input_dir = parent_dir + 'datasets/statistical_dataset/non_functional_dataset/'
phylo_dir = parent_dir + 'datasets/phylogenetic_dataset/'

def plot_distribution(df, width, color, figsize, title):

    plt.style.use('seaborn-white')
    # Data
    ppi = df['prey'].tolist()
    anc = np.array(df['ancestral'].tolist())
    homolog = np.array(df['homolog'].tolist())
    non_homolog = np.array(df['non_homolog'])

    y = np.arange(3)

    fig, ax = plt.subplots(ncols=1, figsize=figsize)

    ax.barh(y + width, anc, width, color='black', edgecolor="black")
    ax.barh(y, homolog, width, color='#939598', edgecolor="black")
    ax.barh(y - width, non_homolog, width, color='white', edgecolor="black")
    ax.set_xlim([0, 60])
    ax.set_yticks(y)
    ax.set_yticklabels(ppi, weight='bold')
    # ax.set_xticks(xticks)
    ax.tick_params(axis='x', direction='out', length=8, width=4, colors=color, labelsize=18)
    ax.tick_params(axis='y', direction='out', length=0, width=0, colors=color, labelsize=18)
    ax.set_xlabel('percentage of non-contingent mutations', fontsize=20, color=color, weight='bold')
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

baits = ['MYO3', 'MYO5']


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

    dict_fisher_exact_anc = dict()
    dict_fisher_exact_hom = dict()
    dict_fisher_exact_anc_non = dict()

    for bait in baits:
        anc_p_list = []
        homolog_p_list = []
        non_homolog_p_list = []

        for prey in preys:

            df1 = pd.read_csv(input_dir + prey + ' - ' + bait + '.dat', sep='\t')
            df = pd.read_csv(variant_dir + prey + ' - ' + bait + '.dat', sep='\t')


            df1 = df1[~((df1['sequence_type_3'] == 'synonymous') & (df1['sequence_type_3'] == 'synonymous'))].copy()
            df1 = df1[~((df1['sequence_type_3'] == 'STOP') & (df1['sequence_type_3'] == 'STOP'))].copy()

            anc_df = df[df['mutated_reference'].isin(anc_list)].copy()
            homolog_df = df[df['mutated_reference'].isin(homolog_list)].copy()
            non_homolog_df = df[~df['mutated_reference'].isin(all_homolog_list)].copy()
            non_homolog_df = non_homolog_df[~non_homolog_df['mutated_reference'].isin(extant_list)].copy()

            anc_df1 = df1[df1['mutated_reference'].isin(anc_list)].copy()
            homolog_df1 = df1[df1['mutated_reference'].isin(homolog_list)].copy()
            non_homolog_df1 = df1[~df1['mutated_reference'].isin(all_homolog_list)].copy()
            non_homolog_df1 = non_homolog_df1[~non_homolog_df1['mutated_reference'].isin(extant_list)].copy()


            perc_anc = round((len(anc_df1.index) / len(anc_df.index)) * 100, 2)
            perc_homolog = round((len(homolog_df1.index) / len(homolog_df.index)) * 100, 2)
            perc_non_homolog = round((len(non_homolog_df1.index) / len(non_homolog_df.index)) * 100, 2)

            anc_p_list.append(perc_anc)
            homolog_p_list.append(perc_homolog)
            non_homolog_p_list.append(perc_non_homolog)

            anc_val = len(anc_df1.index)
            hom_val = len(homolog_df1.index)
            non_val = len(non_homolog_df1.index)

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

        plot_distribution(df, 0.2, '#3b444b', (6, 5), bait)

        list_fe = [dict_fisher_exact_anc, dict_fisher_exact_anc_non, dict_fisher_exact_hom]
        df_fe = pd.DataFrame(list_fe)




