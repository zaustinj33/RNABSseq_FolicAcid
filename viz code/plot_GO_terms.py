import glob
import json, os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#%%
def filter_GO(file):
    # GO terms with a BH corrected p value < 0.01, enrichment score >1.5 and gene counts >10
    test = pd.read_csv(file)#, sep='\t')
    test['Term'] = test['Term'].replace('^.*~','',regex=True)
    test_filter = test.copy() #[((test['Count'] > 3) & (test['Fold Enrichment']) > .5).copy() &
                        #(test['Category'] == 'GOTERM_BP_DIRECT')].copy()# &
                       #(test['FDR'] < 0.1)].copy()
    test_filter['-logBH'] = -np.log10(test_filter['FDR']) #['-logP'] = -np.log10(test_filter['PValue'])

    return test_filter

#%%

#import all sites
#conditions = ["tLF","tMF","tHF", "pLF","pMF","pHF"]
#conditions = ["LFtotal","MFtotal","HFtotal"]  #Jonah files - all m5C
#conditions = ["LFpoly","MFpoly","HFpoly"]  #Jonah files - all m5C
conditions = ["LFMFdown", "HFMFup", "HFMFdown"]  #Jonah files - DTGs
#conditions = ["pMFpLF", "pMFpHF"]  #Jonah files - poly DMS
#conditions = ['tLF_pLF', 'tMF_pMF', 'tHF_pHF']
#conditions = ["LF","MF","HF", 'equalTotal']#, "pLF","pMF","pHF",'equalPoly']
#conditions = ["LFup","LFdown","HFup","HFdown"]
#conditions = ["LFposCorr","LFDMSonly","LFnegCorr","MFposCorr","MFDMSonly","MFnegCorr","HFposCorr","HFDMSonly","HFnegCorr"]

GO_dict = {}
for set in conditions:
    for file in glob.glob("**/"+set+"*GO_anno.csv", recursive=True):
        print(file)
        GO_dict[set] = filter_GO(file)


plot_GO = pd.concat(GO_dict, axis=0).reset_index()
#%%
fig, (ax1) = plt.subplots(1,figsize=(3,6))
#sns.set(rc={'figure.figsize':(3,10)})
sns.set_style("whitegrid")
sns.scatterplot(data=plot_GO, x='level_0', y='Term', size='Count', hue='-logBH',  #logBH
                sizes=(10,200), palette='coolwarm', ax=ax1)
plt.margins(x=0.25)

#ax2.legend(loc='upper right', bbox_to_anchor=(1.7,1), labelspacing=2,
 #          fontsize=14, frameon=False, markerscale=1)
#ax1.tick_params(axis='x', which='major', pad=150)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
plt.xticks(rotation=45)
plt.savefig("GO_bubbleplot_DTGs.png",bbox_inches='tight', dpi=400, transparent=True)
plt.show()
