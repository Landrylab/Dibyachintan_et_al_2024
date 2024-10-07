import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

input_dir = '../Dibyachintan_et_al_2024/datasets/dms_functional_score/'
figure_dir = '../Dibyachintan_et_al_2024/figures/figure_3/'
rsa_dir = '../Dibyachintan_et_al_2024/evo_rate4site/'

df_rsa = pd.read_csv(rsa_dir + 'Myo3.dat', delim_whitespace=True)

myo3_sh3_seq = 'PKFEAAYDFPGSGSSSELPLKKGDIVFISRDEPSGWSLAKLLDGSKEGWVPTAYMTPYK'
preys = ['Bbc1', 'Osh2', 'Pan1', 'Aim21', 'Arc18']

df_list = []
for prey in preys:
    df1 = pd.read_csv(input_dir + prey + '_MYO3_Myo3.dat', sep='\t')

    df1['mutated_reference'] = df1['position'].astype(str) + df1['mut_aa']
    df1 = df1[df1['sequence_type'] != 'synonymous'].copy()
    df1 = df1[df1['sequence_type'] != 'STOP'].copy()
    df1['dF'] = np.clip(df1['dF'].values, a_min=0.0, a_max=1.3)

    df_list.append(df1)

df = pd.concat(df_list)


list_residue = []
for index, row in df_rsa.iterrows():
    evo = float(row['SCORE'])
    if evo <= -1.0:
        list_residue.append('0')
    elif evo <= 0.0:
        list_residue.append('1')
    else:
        list_residue.append('2')

df_rsa['type'] = list_residue
df_rsa = df_rsa[['POS', 'SEQ', 'SCORE', 'type']].copy()
df_rsa.columns = ['position', 'wt_aa', 'evo_rate', 'type']

df_all = pd.merge(df, df_rsa, on=['position', 'wt_aa'])
df_all = df_all.sort_values(by=['type', 'interaction'], ascending=[True, True]).copy()


fig, ax = plt.subplots(1, 1, figsize=(14, 9))
pal = sns.color_palette('Greys', n_colors=7)
sns.boxplot(data=df_all, x="type", y="dF", hue="interaction",
            ax=ax, palette=pal, showfliers=False, linewidth=4)

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.09),
          fancybox=True, shadow=True, ncol=5, prop={'size': 16})

xticks = ['invariant\n(7)', 'slow\n(28)', 'fast\n(24)']
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(labels=xticks, weight='bold')
ax.set_yticks([0.0, 0.5, 1.0, 1.5])
ax.set_yticklabels([0.0, 0.5, 1.0, 1.5], weight='bold')
ax.tick_params(axis='both', which='major', labelsize=24, length=10, width=4, color='black')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_edgecolor('black')
ax.spines['bottom'].set_edgecolor('black')
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.set_ylim([-0.1, 1.5])
ax.set_ylabel('${\Delta}F$', weight='bold', fontsize=24)
ax.set_xlabel('relative rate of substitution', weight='bold', fontsize=24)
plt.savefig(figure_dir + 'figure_3D' + '.pdf', dpi=600, bbox_inches='tight')
plt.show()

