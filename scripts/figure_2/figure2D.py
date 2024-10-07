import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from numpy import mean
from scipy.stats import ttest_rel

parent_dir = '../Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/dms_variant_classification/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_2/'
stat_dir = '../Dibyachintan_et_al_2024/statistical_tests/figure_2/'

preys = ['Pan1', 'Osh2', 'Bbc1', 'Arc18', 'Aim21']

df_list = []
list_p_val = []
list_t_stat = []
list_std = []
list_mean = []
list_cohend = []

for prey in preys:
    df3 = pd.read_csv(input_dir + prey + '_MYO3_Myo3.dat', sep='\t')
    df5 = pd.read_csv(input_dir + prey + '_MYO5_Myo5.dat', sep='\t')

    df3['mutated_reference'] = df3['position'].astype(str) + df3['mut_aa']
    df5['mutated_reference'] = df5['position'].astype(str) + df5['mut_aa']

    df3 = df3[['sequence_type', 'mutated_reference', 'dF', 'classification']].copy()
    df5 = df5[['sequence_type', 'mutated_reference', 'dF', 'classification']].copy()

    #remove synonymous mutations
    df3 = df3[df3['sequence_type'] != 'synonymous'].copy()
    df5 = df5[df5['sequence_type'] != 'synonymous'].copy()

    #cap the upper an lower values of the scores between the median of nonsense and synonymous mutations since vast majority of mutants lie in this dynamic range
    df3['dF'] = np.clip(df3['dF'].values, a_min=0.0, a_max=None)
    df5['dF'] = np.clip(df5['dF'].values, a_min=0.0, a_max=None)

    df = pd.merge(df3, df5, on='mutated_reference')

    list_ddF = []
    for index, row in df.iterrows():

        if row['classification_x'] == 'non-binding' and row['classification_y'] == 'non-binding':
            ddF = row['dF_y'] - row['dF_x']

        else:
            ddF = row['dF_y'] - row['dF_x']
        list_ddF.append(ddF)

    df['ddF'] = list_ddF
    df1 = pd.DataFrame()

    df1['prey'] = [prey] * len(df.index)
    df1['SH3 domain'] = ['Myo3'] * len(df.index)
    df1['ddF'] = df['ddF']
    df_list.append(df1)

    delta_distribution = np.array(df1['ddF'].tolist())

    pval_list = []
    x = np.array(df['dF_x'].tolist())
    y = np.array(df['dF_y'].tolist())

    test = ttest_rel(x, y, alternative='two-sided')
#   test = ttest_1samp(delta_distribution, 0.0, alternative='two-sided')
    t_stat = test.statistic
    pval= test.pvalue
    pval_list.append(pval)

    stdev = np.std(delta_distribution)
    avg = np.mean(delta_distribution)

    #calculate the effect size for difference in functional effects
    cohen_d = round(abs(avg / stdev), 2)

    list_p_val.append(pval)
    list_std.append(stdev)
    list_mean.append(avg)
    list_cohend.append(cohen_d)



df_all = pd.concat(df_list)

df_pval = pd.DataFrame()
df_pval['prey'] = preys
df_pval['p_val'] = list_p_val
df_pval['std'] = list_std
df_pval['mean'] = list_mean
df_pval['cohen_d'] = list_cohend

df_pval.to_excel(stat_dir + 'figure_2D.xlsx', index=False, header=True)


fig, ax = plt.subplots(1, 1, figsize=(8, 8))
res = sns.violinplot(data=df_all, y="prey", x="ddF", ax=ax,  color="silver", cut=0.0)
sns.pointplot(y="prey", x="ddF", data=df_all, estimator=mean, color="black")
ax.axvline(x=0.0, ls='--', color='black')

yticks = preys
ax.set_yticks([0, 1, 2, 3, 4])
ax.set_yticklabels(labels=yticks, weight='bold')
ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax.set_xticklabels([-1.0, -0.5, 0.0, 0.5, 1.0], weight='bold')
ax.tick_params(axis='both', which='major', labelsize=20, length=10, width=4, color='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_edgecolor('black')
ax.spines['bottom'].set_edgecolor('black')
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.set_xlim([-1.0, 1.0])
ax.spines['left'].set_bounds(0, 4)

ax.spines['bottom'].set_position(('axes', -0.01))
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('axes', -0.01))

ax.set_xlabel('difference in functional effects', weight='bold', fontsize=20)
ax.set_ylabel('', weight='bold', fontsize=12)
ax.set_title('${{\Delta}F_{MYO5}^{Myo5}} - {{\Delta}F_{MYO3}^{Myo3}}$', weight='bold', fontsize=24)
plt.savefig(figure_dir + 'figure2D.pdf', dpi=600, bbox_inches='tight')
plt.show()
