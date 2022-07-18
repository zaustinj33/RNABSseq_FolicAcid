from scipy import stats
import os
import pandas as pd
import numpy as np
#BH adjusted Fisher Exact test
#%%

#import all sites
conditions = ['MF', 'KCl0hBS','P0_WT', 'P17_WT', '6W_WT']

#os.chdir('IdeaProjects/FA_project')
total_df = pd.read_csv("allm5C_libraries_filteredDepthAnno.csv", low_memory=False)  # file
#%%

def write_DMS_input(set1, set2, alpha):
    """

    :param set1: prefix of condition 1 from conditions list
    :param set2: prefix of condition 2
    :param alpha: p-value cutoff of FDR test
    :return: DataFrame of sites labeled as DMS or not, with pval and padj values
    """
    overlap1 = pd.read_csv(set1+'_alloverlap.csv')['group']
    overlap2 = pd.read_csv(set2+'_alloverlap.csv')['group']
    group_list = overlap1.tolist() + overlap2.tolist()

    # Aggregate methylation level for each condition (sum counts, divide by total)
    subset_df = total_df[total_df['group'].isin(group_list)]
    subset_df_copy = subset_df.copy()
    subset_df = subset_df.fillna(0)

    regex_input = set1+'|'+set2
    master_set = subset_df.filter(regex=regex_input)
    master_set.index = subset_df['group']


    # aggregate counts and coverage
    count_set = {set1: master_set['C_count_'+set1+'_rep1'] + master_set['C_count_'+set1+'_rep2'],
                 set2: master_set['C_count_'+set2+'_rep1'] + master_set['C_count_'+set2+'_rep2']}

    cov_set = {set1: master_set['cov_'+set1+'_rep1'] + master_set['cov_'+set1+'_rep2'],
               set2: master_set['cov_'+set2+'_rep1'] + master_set['cov_'+set2+'_rep2']}

    ML_dict = {}
    counter = 0

    for i,j in cov_set.items():

        ML_dict[i] = count_set[i].divide(j, fill_value=0)
    pvals = []
    p_adj = []


    try:
        print(len(count_set[set1]) == len(cov_set[set2]))
    except:
        print('data is not same size')
    for i in range(len(count_set[set1])):
        counter += 1
        cont_table = pd.DataFrame({set1:[count_set[set1].iloc[i], cov_set[set1].iloc[i]],
                                   set2:[count_set[set2].iloc[i], cov_set[set2].iloc[i]]})

        odds, pvalue = stats.fisher_exact(cont_table)
        pvals.append(pvalue)
    print(counter)
    pvals_sorted = sorted(pvals, key=float)     # sorted pvalues
    subset_df_copy['pval'] = pvals
    subset_df_copy = subset_df_copy.sort_values('pval', ascending=True)

    rank=1
    for p in pvals_sorted:
        fdr_pval = p*len(pvals_sorted)/rank
        rank+=1
        p_adj.append(fdr_pval)

    subset_df_copy['BH'] = p_adj
    subset_df_copy['shape'] = np.where(subset_df_copy['BH'] <= alpha, 'sig', 'non-sig')
    print("Number of DMS sites: " + str(len(subset_df_copy['group'][subset_df_copy['shape'] == 'sig'])))
    subset_df_copy.to_csv('Compare_'+ set1 + set2 + '_TrueOLtable005.csv')

    return subset_df_copy

#%%

test = pd.DataFrame(write_DMS_input('P17_WT', '6W_WT', 0.05))


