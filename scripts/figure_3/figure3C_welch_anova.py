import pandas as pd
import numpy as np
import scipy.stats as st
import scikit_posthocs as sp
import pingouin as pg

input_dir ='../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
posthoc_dir = '../Dibyachintan_et_al_2024/statistical_tests/figure_3/figure_3C/figure_3C_games_howell/'
stat_dir = '../Dibyachintan_et_al_2024/statistical_tests/figure_3/figure_3C/'
rsa_dir = '../Dibyachintan_et_al_2024/rsa_dssp/'

df_rsa = pd.read_csv(rsa_dir + 'Myo3_rsa.dat', sep='\t')

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
preys = ['Pan1', 'Aim21', 'Arc18', 'Bbc1', 'Osh2']
columns = ['position', 'wt_aa', 'mutated_reference', 'dF', 'interaction']


df1_list = []
for prey in preys:
    df1 = pd.read_csv(input_dir + prey + '_MYO3_Myo3.dat', sep='\t')

    df1['mutated_reference'] = df1['position'].astype(str) + df1['mut_aa']
    df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
    df1 = df1[df1['sequence_type'] != 'STOP'].copy()


    df1 = df1[columns].copy()
    df1.columns = ['position', 'wt_aa', 'mutated_reference', 'dF', 'interaction']
    df1 = df1[columns].copy()
    df1.columns = ['position', 'wt_aa', 'mutated_reference', 'dF', 'interaction']
    df1_list.append(df1)

df2 = pd.concat(df1_list)

list_residue = []
for index, row in df_rsa.iterrows():
    rsa = row['rsa']
    if rsa <= 0.04:
        list_residue.append('B')
    elif rsa <= 0.25:
        list_residue.append('b')
    elif rsa <= 0.50:
        list_residue.append('e')
    else:
        list_residue.append('E')

df_rsa['type'] = list_residue
df_rsa = df_rsa[['resnum', 'aa', 'rsa', 'type']].copy()
df_rsa.columns = ['position', 'wt_aa', 'rsa', 'rsa_type']
df_rsa = df_rsa[(df_rsa['position'] > 0) & (df_rsa['position'] < 60)].copy()

df_all = pd.merge(df2, df_rsa, on=['position', 'wt_aa'])

list_rsa = []
list_stat = []
list_pval = []

rsa_types = ['B', 'b', 'e', 'E']
rsa_type_names = ['completely_buried', 'partially_buried',
                  'partially_exposed', 'completely_exposed']

x_list = []

for index in range(len(rsa_types)):
    rsa_type = rsa_types[index]
    df_t = df_all[df_all['rsa_type'] == rsa_type].copy()
    x = pg.welch_anova(dv='dF', between='interaction', data=df_t)
    x_list.append(x)

    pval_df = pg.pairwise_gameshowell(data=df_t, dv='dF',
                                      between='interaction').round(2)
    pval_df.to_excel(posthoc_dir + rsa_type_names[index] + '_MYO3_Myo3.xlsx', index=True, header=True)

df_x = pd.concat(x_list)
df = pd.DataFrame()
df['position'] = rsa_type_names
df['stat'] = df_x['F'].tolist()
df['pval'] = df_x['p-unc'].tolist()
pvalues = np.array(df['pval'].tolist())
adj_pvalues = p_adjust_bh(pvalues)
df['pval_BH'] = adj_pvalues

df.to_excel(stat_dir + 'figure3C_Welch_ANOVA.xlsx', index=False, header=True)