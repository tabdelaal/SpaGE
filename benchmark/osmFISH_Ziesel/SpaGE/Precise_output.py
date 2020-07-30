import os
os.chdir('osmFISH_Ziesel/')

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(1,'SpaGE/')
from principal_vectors import PVComputation

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data = datadict['osmFISH_data']
osmFISH_data_scaled = datadict['osmFISH_data_scaled']
osmFISH_meta= datadict['osmFISH_meta']
del datadict

with open ('data/SpaGE_pkl/Ziesel.pkl', 'rb') as f:
    datadict = pickle.load(f)

RNA_data = datadict['RNA_data']
RNA_data_scaled = datadict['RNA_data_scaled']
RNA_meta = datadict['RNA_meta']
del datadict

Common_data = RNA_data_scaled[np.intersect1d(osmFISH_data_scaled.columns,RNA_data_scaled.columns)]

n_factors = 30
n_pv = 30
n_pv_display = 30
dim_reduction = 'pca'
dim_reduction_target = 'pca'

pv_FISH_RNA = PVComputation(
        n_factors = n_factors,
        n_pv = n_pv,
        dim_reduction = dim_reduction,
        dim_reduction_target = dim_reduction_target
)

pv_FISH_RNA.fit(Common_data,osmFISH_data_scaled[Common_data.columns])

fig = plt.figure()
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('osmFISH',fontsize=18, color='black')
plt.ylabel('Ziesel_RNA',fontsize=18, color='black')
plt.xticks(np.arange(n_pv_display)+0.5, range(1, n_pv_display+1), fontsize=12)
plt.yticks(np.arange(n_pv_display)+0.5, range(1, n_pv_display+1), fontsize=12, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.show()

plt.figure()
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=14,color='black')
    
plt.xlabel('osmFISH',fontsize=18, color='black')
plt.ylabel('Ziesel_RNA',fontsize=18, color='black')
plt.xticks(np.arange(n_pv_display)+0.5, range(1, n_pv_display+1), fontsize=12)
plt.yticks(np.arange(n_pv_display)+0.5, range(1, n_pv_display+1), fontsize=12, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:30]

### Technology specific Processes
Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)

# explained variance RNA
np.sum(pv_FISH_RNA.source_explained_variance_ratio_[np.arange(Effective_n_pv)])*100
# explained variance spatial
np.sum(pv_FISH_RNA.target_explained_variance_ratio_[np.arange(Effective_n_pv)])*100
