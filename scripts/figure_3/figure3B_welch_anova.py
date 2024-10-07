import pandas as pd
import numpy as np
import scipy.stats as st
import scikit_posthocs as sp
import pingouin as pg

input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
posthoc_dir = '../Dibyachintan_et_al_2024/statistical_tests/figure_3/figure_3B/figure_3B_games_howell/'
stat_dir = '../Dibyachintan_et_al_2024/statistical_tests/figure_3/figure_3B/'

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

df_list = []

df = pd.DataFrame()
columns = ['position', 'mutated_reference', 'dF', 'interaction']
count = 0
df1_list = []
for prey in preys:
    df1 = pd.read_csv(input_dir + prey + '_MYO3_Myo3.dat', sep='\t')

    df1['mutated_reference'] = df1['position'].astype(str) + df1['mut_aa']
    df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
    df1 = df1[df1['sequence_type'] != 'STOP'].copy()
    # df1['dF'] = np.clip(df1['dF'].values, a_min=0.0, a_max=None)

    df1 = df1[columns].copy()
    df1.columns = ['position', 'mutated_reference', 'dF', 'interaction']
    df1_list.append(df1)


df_all = pd.concat(df1_list)

list_pos = []
list_stat = []
list_pval = []
list_pval_bh = []

x_list = []

for pos in range(1, 60):
    df_t = df_all[df_all['position'] == pos].copy()
    x = pg.welch_anova(dv='dF', between='interaction', data=df_t)
    x_list.append(x)
    list_pos.append(pos)


    pval_df = pg.pairwise_gameshowell(data=df_t, dv='dF',
                                      between='interaction').round(3)
    pval_df.to_excel(posthoc_dir + str(pos) + '_MYO3_Myo3.xlsx', index=True, header=True)


df_x = pd.concat(x_list)
df = pd.DataFrame()
df['position'] = list_pos
df['stat'] = df_x['F'].tolist()
df['pval'] = df_x['p-unc'].tolist()
pvalues = np.array(df['pval'].tolist())
adj_pvalues = p_adjust_bh(pvalues)
df['pval_BH'] = adj_pvalues
df2 = df[df['pval_BH'] <= 0.01]
df.to_excel(stat_dir + 'figure3B_Welch_ANOVA.xlsx', index=False, header=True)