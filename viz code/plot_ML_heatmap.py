import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
from sklearn.cluster import AgglomerativeClustering

#%%
#import all sites
conditions = ["tLF","tMF","tHF","pLF","pMF","pHF"]

total_file = "allm5C_libraries_filteredDepthAnno.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file

# Subset by gene name if needed
gene_list = ['RPS26', 'RPL41', 'RPL18A', 'RPLP1', 'RPL36', 'RPL18'
             ]
#total_df = total_df[total_df['gene_name'].str.upper().isin(gene_list)]
#%%

# Aggregate methylation level for each condition (sum counts, divide by total)
total_df.index = total_df['group']
cov_df = total_df.filter(regex='cov')
count_df = total_df.filter(regex='count')

cov_dict = {}
count_dict = {}
for name in conditions:
    cov_dict[name] = cov_df.filter(regex=name).sum(axis=1)
    count_dict[name] = count_df.filter(regex=name).sum(axis=1)

ML_dict = {}

for i,j in cov_dict.items():
    ML_dict[i] = count_dict[i].divide(j, fill_value=0)


result_df = pd.DataFrame(ML_dict)

result_df.replace(np.nan, 0, inplace=True)
#result_df.replace(0, np.nan, inplace=True)

result_df = result_df[(result_df['pLF'] >= 0.1) | (result_df['pMF'] >= 0.1) | (result_df['pHF'] >= 0.1) |
                      (result_df['tLF'] >= 0.1) | (result_df['tMF'] >= 0.1) | (result_df['tHF'] >= 0.1)]

result_df.dropna(axis=0, inplace=True)
test = result_df[result_df.index.isin(total_df['group'])]
result_df.to_csv("AllOverlap_methylationLevel_total.csv")

#%%
from matplotlib.colors import LinearSegmentedColormap

boundaries = [0.0, 0.05, 0.1, 0.2, 0.4 ,0.6,1.0]
hex_colors = sns.color_palette("RdYlBu_r", n_colors=len(boundaries) * 2).as_hex()
hex_colors = [hex_colors[i] for i in range(0, len(hex_colors), 2)]

colors=list(zip(boundaries, hex_colors))

custom_color_map = LinearSegmentedColormap.from_list(
    name="cus",
    colors=colors,
)
# %%
# Define clusters
def generate_clusters(input_df, n_clusters, method):
    correlations_array = np.asarray(input_df)

    row_linkage = hierarchy.linkage(
        distance.pdist(correlations_array), method=method)

    col_linkage = hierarchy.linkage(
        distance.pdist(correlations_array.T), method=method)

    model = AgglomerativeClustering(n_clusters=n_clusters, affinity='euclidean', linkage=method)
    model = model.fit_predict(correlations_array)

    return model, row_linkage, col_linkage


# %%
n_clusters = 6
model_list = generate_clusters(result_df, n_clusters, 'ward')

# %%
cluster_palette = ['red','blue','green','orange','purple','pink']
lut = dict(zip(set(model_list[0]), cluster_palette))
row_colors = pd.DataFrame(model_list[0])[0].map(lut)

cg=sns.clustermap(result_df.reset_index(drop=True), row_linkage=model_list[1], col_linkage=model_list[2], cmap=custom_color_map,
                  row_colors = row_colors, figsize=(5,5), yticklabels=False, col_cluster=False,
                  robust=True, method='ward')  #, row_cluster=False)  # z_score=0,
cg.ax_row_dendrogram.set_visible(False)
#plt.savefig("ML_conditions_clusteringHeatmapDepthTotal.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()

#%%
# Annotate clusters
result_df['clusterID'] = model_list[0]
total_df.index = total_df['group']

All_cluster = total_df.merge(result_df,left_index=True, right_index=True, how='left')
All_cluster.dropna(axis = 0, how = 'any', inplace = True, subset=['clusterID'])

All_cluster.to_csv("AllOverlap_methylationLevel_clustersTotal.csv")
#%%
from scipy.stats import zscore
# write correlation matrix (z-score)

zscore_vals = result_df.apply(zscore, axis=1)

#%%
# Correlation matix of samples
from scipy.spatial import distance
from scipy.cluster import hierarchy

correlations = result_df.corr()
correlations_array = np.asarray(result_df.corr())

row_linkage = hierarchy.linkage(
    distance.pdist(correlations_array), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(correlations_array.T), method='average')

sns.clustermap(correlations, row_linkage=col_linkage, col_linkage=row_linkage, method="average",
               figsize=(5, 10))
plt.show()


#%%

# save cluster information

total_clusters = pd.read_csv('AllOverlap_methylationLevel_clusters_Poly.csv')
def save_cluster(name, clusterID):
    cluster = total_clusters[total_clusters['clusterID'].isin(clusterID)]
    cluster_IDs = cluster['gene_name'].unique()
    np.savetxt(name+"_cluster_IDs.txt", cluster_IDs, fmt='%s')

#%%
save_cluster('equalPoly', [6,3,5])

#%%
# DMS heatmap
## 4 Plot
plot_LF = pd.read_csv('LFC_ML_compare_LF005.csv')
plot_MF = pd.read_csv('LFC_ML_compare_MF005.csv')
plot_HF = pd.read_csv('LFC_ML_compare_HF005.csv')
#%%
plot_polyDMS_LF = pd.read_csv('Compare_pMFpLF_TrueOLtable005.csv')
plot_polyDMS_LF['ML_1'] = (plot_polyDMS_LF['C_count_pMF_rep1'] + plot_polyDMS_LF['C_count_pMF_rep2'])/ \
                          (plot_polyDMS_LF['cov_pMF_rep1'] + plot_polyDMS_LF['cov_pMF_rep2'])
plot_polyDMS_LF['ML_2'] = (plot_polyDMS_LF['C_count_pLF_rep1'] + plot_polyDMS_LF['C_count_pLF_rep2'])/\
                          (plot_polyDMS_LF['cov_pLF_rep1'] + plot_polyDMS_LF['cov_pLF_rep2'])

plot_polyDMS_HF = pd.read_csv('Compare_pMFpHF_TrueOLtable005.csv')
plot_polyDMS_HF['ML_1'] = (plot_polyDMS_HF['C_count_pMF_rep1'] + plot_polyDMS_HF['C_count_pMF_rep2'])/ \
                          (plot_polyDMS_HF['cov_pMF_rep1'] + plot_polyDMS_HF['cov_pMF_rep2'])
plot_polyDMS_HF['ML_2'] = (plot_polyDMS_HF['C_count_pHF_rep1'] + plot_polyDMS_HF['C_count_pHF_rep2'])/ \
                          (plot_polyDMS_HF['cov_pHF_rep1'] + plot_polyDMS_HF['cov_pHF_rep2'])
#%%

# Define clusters
def subset_df_to_resultsDF(DMS_DF, name, hypo_hyper):
    results_DMS_all = DMS_DF[DMS_DF['shape'] == 'sig']
    results_DMS_df = results_DMS_all[['gene_name', 'methRate_pMF_rep1', 'methRate_pMF_rep2',
                                      'methRate_'+name+'_rep1','methRate_'+name+'_rep2']]
    results_DMS_df = results_DMS_df.set_index('gene_name')
    #results_DMS_df = results_DMS_df.dropna()

    if hypo_hyper == 'hypo':
        results_DMS_df = results_DMS_df[results_DMS_df['methRate_pMF_rep1'] > results_DMS_df['methRate_'+name+'_rep1']]
        results_DMS_df = results_DMS_df.sort_values('methRate_pMF_rep1', ascending=False)
    else:
        results_DMS_df = results_DMS_df[results_DMS_df['methRate_pMF_rep1'] < results_DMS_df['methRate_'+name+'_rep1']]
        results_DMS_df = results_DMS_df.sort_values('methRate_'+name+'_rep1', ascending=False)
    return results_DMS_df

plot_DMS_df = subset_df_to_resultsDF(plot_polyDMS_HF, 'pHF', 'hyper')

# if needed
#model_list = generate_clusters(plot_DMS_df, 2, 'single')

#%%
cg=sns.clustermap(plot_DMS_df,  #row_linkage=row_linkage,
                  figsize=(5,5), col_cluster=False, row_cluster=False, cmap=custom_color_map,
                  robust=True, method='single')  #, row_cluster=False)  # z_score=0,
cg.ax_row_dendrogram.set_visible(False)
plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0)
#plt.savefig("figures/Heatmap_DMS_PolyHF_hyper005.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()


#%%
# Total v Poly heatmaps
plot_polyDMS_LF = pd.read_csv('Compare_tLFpLF_TrueOLtable005.csv')
plot_polyDMS_LF['ML_1'] = (plot_polyDMS_LF['C_count_pLF_rep1'] + plot_polyDMS_LF['C_count_pLF_rep2'])/ \
                          (plot_polyDMS_LF['cov_pLF_rep1'] + plot_polyDMS_LF['cov_pLF_rep2'])
plot_polyDMS_LF['ML_2'] = (plot_polyDMS_LF['C_count_tLF_rep1'] + plot_polyDMS_LF['C_count_tLF_rep2'])/ \
                          (plot_polyDMS_LF['cov_tLF_rep1'] + plot_polyDMS_LF['cov_tLF_rep2'])

plot_polyDMS_HF = pd.read_csv('Compare_pMFpHF_TrueOLtable005.csv')
plot_polyDMS_HF['ML_1'] = (plot_polyDMS_HF['C_count_pMF_rep1'] + plot_polyDMS_HF['C_count_pMF_rep2'])/ \
                          (plot_polyDMS_HF['cov_pMF_rep1'] + plot_polyDMS_HF['cov_pMF_rep2'])
plot_polyDMS_HF['ML_2'] = (plot_polyDMS_HF['C_count_pHF_rep1'] + plot_polyDMS_HF['C_count_pHF_rep2'])/ \
                          (plot_polyDMS_HF['cov_pHF_rep1'] + plot_polyDMS_HF['cov_pHF_rep2'])
#%%

# Define clusters
def subset_df_to_resultsDF(DMS_DF, name, hypo_hyper):
    results_DMS_all = DMS_DF[DMS_DF['shape'] == 'sig']
    results_DMS_df = results_DMS_all[['gene_name', 'methRate_pMF_rep1', 'methRate_pMF_rep2',
                                      'methRate_'+name+'_rep1','methRate_'+name+'_rep2']]
    results_DMS_df = results_DMS_df.set_index('gene_name')
    #results_DMS_df = results_DMS_df.dropna()

    if hypo_hyper == 'hypo':
        results_DMS_df = results_DMS_df[results_DMS_df['methRate_pMF_rep1'] > results_DMS_df['methRate_'+name+'_rep1']]
        results_DMS_df = results_DMS_df.sort_values('methRate_pMF_rep1', ascending=False).head(n=3)

    else:
        results_DMS_df = results_DMS_df[results_DMS_df['methRate_pMF_rep1'] < results_DMS_df['methRate_'+name+'_rep1']]
        results_DMS_df = results_DMS_df.sort_values('methRate_'+name+'_rep1', ascending=False).head(n=7)
    return results_DMS_df

plot_DMS_df = subset_df_to_resultsDF(plot_polyDMS_HF, 'pHF', 'hyper')

# if needed, clustering
#model_list = generate_clusters(plot_DMS_df, 2, 'single')

#%%
cg=sns.clustermap(plot_DMS_df,  #row_linkage=row_linkage,
                  figsize=(5,5), col_cluster=False, row_cluster=False, cmap=custom_color_map,
                  robust=True, method='single')#, yticklabels=False)  # z_score=0,
cg.ax_row_dendrogram.set_visible(False)
plt.setp(cg.ax_heatmap.get_yticklabels(), rotation=0)
plt.savefig("figures/Heatmap_DMS_HF_hyper005.png", bbox_inches='tight', dpi=400, transparent=True)
plt.show()

plt.close()