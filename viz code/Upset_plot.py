import os, glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from upsetplot import UpSet

#%%
# 1) import
conditions = ["tLF","tMF","tHF", "pLF","pMF","pHF"]

m5Csite_dict = {}
for file in glob.glob("**/*_allOverlap.csv", recursive=True):
    name = file.strip('_allOverlap.csv')
    m5Csite_dict[name] = pd.read_csv(file)['group']
    print(name)

#%%
from upsetplot import from_contents
m5Cplot_format = from_contents(m5Csite_dict)
m5Cplot_format.head()

#%%
p = UpSet(m5Cplot_format, subset_size='count', min_subset_size = 3).plot()
#p.add_catplot(value='methRate', kind='strip', color='black')
plt.savefig("OverlapConditions_upsetPlot.png", bbox_inches='tight', dpi=200, transparent=True)

plt.show()


