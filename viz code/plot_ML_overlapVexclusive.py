import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#%%
# Create ML of sites overlapped ([# replicates in set] + 1)
def calcML(df):
    ML_df = df[df.filter(like="methRate").columns].mean(axis=1)
    return ML_df


conditions = ["tLF","tMF","tHF"] #["pLF","pMF","pHF"]#

#%%
#import all sites
total_file = "allm5C_libraries_6Cfiltered.csv"
total_df = pd.read_csv(total_file, low_memory=False)  # file

#MT
total_df['chrom'] = total_df['group'].str.replace(" .*$","")
#total_df = total_df[total_df['chrom'] !=  'MT']

overlapped_site_ML = {}
exclusive_site_ML = {}
for set in conditions:
    print(set)
    ext = "_allOverlap.csv"
    #ext = "_allUnion.csv"
    for file in glob.glob("**/"+set+ext, recursive=True):
        print(file)
        df_set = pd.read_csv(file, low_memory=False)
        name = "_"+set
        df = df_set[df_set.filter(like=name).columns]
        df['group'] = df_set['group']
        df = df[df['group'].isin(total_df['group'])]

        # subset total site columns df by name
        df_exclusive = total_df[total_df.filter(like=name).columns]
        df_exclusive['group'] = total_df['group']
        # Remove sites
        df_exclusive = df_exclusive[~df_exclusive['group'].isin(df_set['group'])]

        print("subset overlapped")
        overlapped_site_ML[set] = calcML(df)
        print("subset exlusive")
        exclusive_site_ML[set] = calcML(df_exclusive)

overlapped_df = pd.DataFrame(overlapped_site_ML)
exclusive_df = pd.DataFrame(exclusive_site_ML)
#%%
flierprops = dict(markersize=5, marker='.',
                  linestyle='none')
palette = ["#00FF00","#0000FF","#FF0000"]
sns.set(rc = {'figure.figsize':(3,4)})
custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)

p = sns.boxplot(data=overlapped_df, linewidth=2, palette=palette, flierprops = flierprops)
#sns.scatterplot(data=overlapped_df, color='black')
plt.ylim(0,1)
plt.savefig("overlapped_ML_boxplot_trueOverlap_total.png", bbox_inches='tight', dpi=200, transparent=True)
plt.show()
#%%
flierprops = dict(markersize=5, marker='.',
                  linestyle='none')
sns.boxplot(data=exclusive_df, linewidth=2,palette='bright', flierprops = flierprops)
plt.ylim(0,1)
plt.savefig("exclusive_ML_boxplot_poly.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()

#%%

from scipy.stats import ranksums
#%%
ranksums(overlapped_df['tMF'].dropna(), overlapped_df['tLF'].dropna())
