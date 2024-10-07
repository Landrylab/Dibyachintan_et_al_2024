import pandas as pd

import matplotlib.pyplot as plt


ax_color = '#037ab1'
fontsize = 24
# ax_color = '#49997c'

parent_dir =  '../Dibyachintan_et_al_2024/'
input_dir = parent_dir + 'datasets/functional_scores_chimera/'
figure_dir = parent_dir + 'figures/figure_6/'

df1 = pd.read_csv(input_dir + 'MYO3_OSH2.dat', sep='\t')
x = 26
y = 39
z = 56

df2 = df1[df1['type'] == 'double_chimera'].copy()
df2 = df2[['type', 'position_2', 'AA_2',
           'position_1', 'AA_1',
           'F1', 'F2', 'F3',
           'F', 'SE_F']].copy()

df2.columns = ['type', 'position_1', 'AA_1',
               'position_2', 'AA_2',
               'F1', 'F2', 'F3',
               'F', 'SE_F']

df_list = [df1, df2]
df = pd.concat(df_list)
df = df.sort_values(by=['type', 'position_1', 'position_2'],
                    ascending=[False, True, True])

myo3_anc = [15, 25, 27, 29, 31, 58]
myo5_anc = [2, 26, 39, 56, 59]

df_wt = df[(df['position_1'] == 0) & (df['position_2'] == 0)].copy()
df_m = df[df['position_1'] == x].copy()
df_nm = df[(df['position_1'] != x) & (df['position_2'] == 0)].copy()

mut_wt = df_wt['F'].tolist()[0]
mut_del = df_m[df_m['position_2'] == 0]['F'].tolist()[0]
mut_y = df_nm[df_nm['position_1'] == y]['F'].tolist()[0]
mut_double_y = df_m[df_m['position_2'] == y]['F'].tolist()[0]

mut_wt_err = df_wt['SE_F'].tolist()[0]
mut_del_err = df_m[df_m['position_2'] == 0]['SE_F'].tolist()[0]
mut_y_err = df_nm[df_nm['position_1'] == y]['SE_F'].tolist()[0]
mut_double_y_err = df_m[df_m['position_2'] == y]['SE_F'].tolist()[0]


mut_wt = df_wt['F'].tolist()[0]
mut_del = df_m[df_m['position_2'] == 0]['F'].tolist()[0]
mut_z = df_nm[df_nm['position_1'] == z]['F'].tolist()[0]
mut_double_z = df_m[df_m['position_2'] == z]['F'].tolist()[0]

mut_wt_err = df_wt['SE_F'].tolist()[0]
mut_del_err = df_m[df_m['position_2'] == 0]['SE_F'].tolist()[0]
mut_z_err = df_nm[df_nm['position_1'] == z]['SE_F'].tolist()[0]
mut_double_z_err = df_m[df_m['position_2'] == z]['SE_F'].tolist()[0]

x_arr = [0, 1, 1, 2]
y_arr = [mut_wt, mut_del, mut_y, mut_double_y]
y_arr_err = [mut_wt_err, mut_del_err, mut_y_err, mut_double_y_err]

x_arr = [0, 1, 1, 2]
z_arr = [mut_wt, mut_del, mut_z, mut_double_z]
z_arr_err = [mut_wt_err, mut_del_err, mut_z_err, mut_double_z_err]


fig, ax = plt.subplots(figsize=(8, 8))
ax.errorbar(x_arr, y_arr, yerr=y_arr_err, fmt='o', capsize=3.0,
            mfc='#231f20', ecolor='#231f20',
            mec='#231f20', ms=5, mew=4)

ax.plot((x_arr[0], x_arr[1]), (y_arr[0], y_arr[1]), 'k-')
ax.plot((x_arr[0], x_arr[2]), (y_arr[0], y_arr[2]), 'k-')

ax.plot((x_arr[1], x_arr[3]), (y_arr[1], y_arr[3]), 'k-')
ax.plot((x_arr[2], x_arr[3]), (y_arr[2], y_arr[3]), 'k-')


ax.errorbar(x_arr, z_arr, yerr=z_arr_err, fmt='o', capsize=3.0,
            mfc='#231f20', ecolor='#231f20',
            mec='#231f20', ms=5, mew=4)

ax.plot((x_arr[0], x_arr[1]), (z_arr[0], z_arr[1]), 'k-')
ax.plot((x_arr[0], x_arr[2]), (z_arr[0], z_arr[2]), 'k-')

ax.plot((x_arr[1], x_arr[3]), (z_arr[1], z_arr[3]), 'k-')
ax.plot((x_arr[2], x_arr[3]), (z_arr[2], z_arr[3]), 'k-')


ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_edgecolor(ax_color)
ax.spines['bottom'].set_edgecolor(ax_color)
ax.spines['left'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)

ax.set_xticks([0, 1, 2])
ax.set_xticklabels(labels=['0', '1', '2'], weight='bold')

# for Osh2
ax.set_yticks([2, 0, -2, -4])
ax.set_yticklabels(labels=['2', '0', '-2', '-4'], weight='bold')


ax.tick_params(axis='both', which='major', labelsize=fontsize, length=10, width=4, color=ax_color)
plt.savefig(figure_dir + '6D.pdf', dpi=600, bbox_inches='tight')
plt.show()




