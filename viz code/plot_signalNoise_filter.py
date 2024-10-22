import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

os.chdir('IdeaProjects/FA_project')
"""
Input files to filter signal/noise and FDR sites:
1) annotated meRanGh output file  (*annotateALL.txt)
2) C-cutoff + standard filter meRanGh file (*3_cutoff_basicFilter.txt)
3) FDR annotated meRanGh output (*FDR_FDR.txt)
"""

#%%
pd.options.mode.chained_assignment = None  # default='warn'
def filter_m5C(input_file):
    df = pd.read_csv(input_file, sep='\t', low_memory=False)  # file
    #1) basic filter (remove (df['cov'] >= 20) & for MT libraries
    df_filter = df[( (df['methRate'] >= 0.1) & (df['C_count'] >= 6) & (df['cov'] >= 20) &
                    (df['state'] == 'M'))]

    # Important columns
    #imp_cols = ['#SeqID', 'cov', 'C_count', 'methRate', 'site_loc']

    # Add unique site info, subset by important columns
    df_filter['site_loc'] = df_filter['#SeqID'].str.replace('chr','') +\
                            "_" + df_filter['refPos'].astype(str)
    df_filter['site_loc'] = df_filter['site_loc'].str.replace('M_','MT_')
    #df_final = df_filter[df_filter.columns.intersection(imp_cols)]

    return df_filter

sample_list = ['LF_rep1', 'LF_rep2', 'MF_rep1', 'MF_rep2', 'HFold_rep1', 'HFold_rep2',
               'pLF_rep1', 'pLF_rep2', 'pMF_rep1', 'pMF_rep2',
                'pHF_rep1', 'pHF_rep2']

#%%

SN_dict = {}
drop_list = {}
for name in sample_list:
    try:
        signal_file = glob.glob("**/"+name+"_Genome10XCall_3_Cutoff_basicFilter.txt", recursive=True)
        noise_file = glob.glob("**/"+name+"_Genome10XCall_annotate.txt", recursive=True)

        noise = filter_m5C(noise_file[0])
        print(noise_file)
        print(len(noise))

        # load files
        signal = filter_m5C(signal_file[0])
        #signal_gene = signal[(signal['gene'] != 'gene:no_feature;')]  # only genes with features
        print(signal_file)
        print(len(signal))


        #1) filter if N(signal)/N(noise) >=0.9; where N is number of reads at a given site
        # .x is signal, .y is noise
        intersect = signal.merge(noise, on='site_loc')
        #print(len(intersect))

        SN = (intersect['cov_x'] + intersect['C_count_x']) / (intersect['cov_y'] + intersect['C_count_y'])
        SN_dict[name] = SN
        noisy_sites = intersect['site_loc'][((intersect['cov_x'] + intersect['C_count_x']) /
                                             (intersect['cov_y'] + intersect['C_count_y']) < 0.9)]
        drop_list[name] = noisy_sites
        print("S/N sites removed:" + str(len(noisy_sites)))
        final_df = signal[~signal['site_loc'].isin(noisy_sites)]
        final_df.to_csv(name+"_Genome10xCall_3signalFilter.txt", index=False, sep='\t')

    except (ValueError, IndexError):
        print(name + " signal and or noise file not found")
        final_df = {}
        pass

    try:

        #2) retain site if p-value adjusted (FDR) is <= 0.05.
        FDR_file = glob.glob("**/"+name+"_Genome10xCall_3_Cutoff_FDR_FDR_0.05.txt", recursive=True)
        FDR = filter_m5C(FDR_file[0])
        intersect_FDR = FDR.merge(final_df, on='site_loc')
        print(len(final_df))

        goodFDR_sites = intersect_FDR['site_loc'][intersect_FDR['p-value_mState_adj'] <= 0.05]
        print("FDR sites removed:" + str(len(final_df) - len(goodFDR_sites)))
        print(len(goodFDR_sites))
        drop_list[name] = goodFDR_sites
        final_df_fdr = final_df[final_df['site_loc'].isin(goodFDR_sites)]
        final_df_fdr.to_csv(name+"_Genome10xCall_3fdrFilter.txt", index=False, sep='\t')

    except (ValueError, IndexError):
        print(name + " FDR file not found")
        pass


#%%
plot_SN_df = pd.DataFrame(SN_dict)
plot_SN_df = plot_SN_df[plot_SN_df <= 1]  # remove outlier points (signal > noise)
#plot_SN_df['group'] = plot_SN_df['name'].str.replace(r'817.*$|_.*$', '')  # add levels

palette_bright = sns.color_palette("bright",12)
color_order = [1,2,3,4,6,7,9,10,11,0,5,8]
colors = [palette_bright[i] for i in color_order]

graph = sns.boxplot(data=plot_SN_df, palette=colors)
graph.axhline(0.9, ls='--', color='black')
plt.savefig("SignalNoise_boxplot.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

binary_SN_DF_fail = plot_SN_df[plot_SN_df < 0.9].count()
binary_SN_DF_pass = plot_SN_df[plot_SN_df >= 0.9].count()

binary_SN_DF_sum = binary_SN_DF_pass + binary_SN_DF_fail
binary_SN_DF_pass = binary_SN_DF_pass / binary_SN_DF_sum
binary_SN_DF_fail = binary_SN_DF_fail / binary_SN_DF_sum

binary_SN_DF = pd.concat([binary_SN_DF_pass, binary_SN_DF_fail], axis=1).reset_index()

binary_SN_DF = binary_SN_DF.rename({'index':'library', 0:'pass', 1:'fail'}, axis=1)
binary_SN_DF = binary_SN_DF[~binary_SN_DF["library"].str.contains("Gall")]

binary_SN_DF_long = pd.melt(binary_SN_DF, id_vars='library', value_vars=['pass','fail'])
#%%

binary_SN_DF.plot(kind='bar', stacked=True, color=['red','blue'],
              linewidth=1, edgecolor='black')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.savefig("SignalNoise_reads_removed.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

