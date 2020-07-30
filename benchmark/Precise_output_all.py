#%% STARmap_AllenVISp
import os
os.chdir('STARmap_AllenVISp/')

import pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(1,'SpaGE/')
from principal_vectors import PVComputation

with open ('data/SpaGE_pkl/Starmap.pkl', 'rb') as f:
    datadict = pickle.load(f)

Starmap_data_scaled = datadict['Starmap_data_scaled']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data_scaled = datadict['RNA_data_scaled']
del datadict

Common_data = RNA_data_scaled[np.intersect1d(Starmap_data_scaled.columns,RNA_data_scaled.columns)]

n_factors = 50
n_pv = 50
n_pv_display = 50
dim_reduction = 'pca'
dim_reduction_target = 'pca'

pv_FISH_RNA = PVComputation(
        n_factors = n_factors,
        n_pv = n_pv,
        dim_reduction = dim_reduction,
        dim_reduction_target = dim_reduction_target
)

pv_FISH_RNA.fit(Common_data,Starmap_data_scaled[Common_data.columns])

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('Starmap',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('Starmap',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:50]

Starmap_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
Starmap_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%% osmFISH_Ziesel
os.chdir('osmFISH_Ziesel/')

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data_scaled = datadict['osmFISH_data_scaled']
del datadict

with open ('data/SpaGE_pkl/Ziesel.pkl', 'rb') as f:
    datadict = pickle.load(f)

RNA_data_scaled = datadict['RNA_data_scaled']
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

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Zeisel',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Zeisel',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:30]

Ziesel_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
Ziesel_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%% osmFISH_AllenVISp
os.chdir('osmFISH_AllenVISp/')

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data_scaled = datadict['osmFISH_data_scaled']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data_scaled = datadict['RNA_data_scaled']
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

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:30]

AllenVISp_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
AllenVISp_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%% osmFISH_AllenSSp
os.chdir('osmFISH_AllenSSp/')

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data_scaled = datadict['osmFISH_data_scaled']
del datadict

with open ('data/SpaGE_pkl/Allen_SSp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data_scaled = datadict['RNA_data_scaled']
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

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Allen_SSp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('osmFISH',fontsize=12, color='black')
plt.ylabel('Allen_SSp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:30]

AllenSSp_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
AllenSSp_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%% MERFISH_Moffit
os.chdir('MERFISH_Moffit/')

with open ('data/SpaGE_pkl/MERFISH.pkl', 'rb') as f:
    datadict = pickle.load(f)

MERFISH_data_scaled = datadict['MERFISH_data_scaled']
del datadict

with open ('data/SpaGE_pkl/Moffit_RNA.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data_scaled = datadict['RNA_data_scaled']
del datadict

Common_data = RNA_data_scaled[np.intersect1d(MERFISH_data_scaled.columns,RNA_data_scaled.columns)]

n_factors = 50
n_pv = 50
n_pv_display = 50
dim_reduction = 'pca'
dim_reduction_target = 'pca'

pv_FISH_RNA = PVComputation(
        n_factors = n_factors,
        n_pv = n_pv,
        dim_reduction = dim_reduction,
        dim_reduction_target = dim_reduction_target
)

pv_FISH_RNA.fit(Common_data,MERFISH_data_scaled[Common_data.columns])

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('MERFISH',fontsize=12, color='black')
plt.ylabel('Moffit',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('MERFISH',fontsize=12, color='black')
plt.ylabel('Moffit',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:50]

MERFISH_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
MERFISH_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%%seqFISH_AllenVISp
os.chdir('seqFISH_AllenVISp/')

with open ('data/SpaGE_pkl/seqFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

seqFISH_data = datadict['seqFISH_data']
seqFISH_data_scaled = datadict['seqFISH_data_scaled']
seqFISH_meta= datadict['seqFISH_meta']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data = datadict['RNA_data']
RNA_data_scaled = datadict['RNA_data_scaled']
del datadict

all_centroids  = seqFISH_meta[['X','Y']]
cmap = cm.get_cmap('viridis',20)

Common_data = RNA_data_scaled[np.intersect1d(seqFISH_data_scaled.columns,RNA_data_scaled.columns)]

n_factors = 50
n_pv = 50
n_pv_display = 50
dim_reduction = 'pca'
dim_reduction_target = 'pca'

pv_FISH_RNA = PVComputation(
        n_factors = n_factors,
        n_pv = n_pv,
        dim_reduction = dim_reduction,
        dim_reduction_target = dim_reduction_target
)

pv_FISH_RNA.fit(Common_data,seqFISH_data_scaled[Common_data.columns])

plt.figure(figsize=(4,4))
sns.heatmap(pv_FISH_RNA.initial_cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
plt.xlabel('seqFISH+',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

plt.figure(figsize=(4, 4))
sns.heatmap(pv_FISH_RNA.cosine_similarity_matrix_[:n_pv_display,:n_pv_display], cmap='seismic_r',
            center=0, vmax=1., vmin=0)
for i in range(n_pv_display-1):
    plt.text(i+1,i+.7,'%1.2f'%pv_FISH_RNA.cosine_similarity_matrix_[i,i], fontsize=5,color='black')
    
plt.xlabel('seqFISH+',fontsize=12, color='black')
plt.ylabel('Allen_VISp',fontsize=12, color='black')
plt.xticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8)
plt.yticks(np.arange(4,n_pv_display,5)+0.5, range(5, n_pv_display+1,5), fontsize=8, rotation='horizontal')
plt.gca().set_ylim([n_pv_display,0])
plt.axis('scaled')
plt.show()

Importance = pd.Series(np.sum(pv_FISH_RNA.source_components_**2,axis=0),index=Common_data.columns)
Importance.sort_values(ascending=False,inplace=True)
Importance.index[0:50]

seqFISH_before = np.diag(pv_FISH_RNA.initial_cosine_similarity_matrix_)
seqFISH_after = np.diag(pv_FISH_RNA.cosine_similarity_matrix_)

#%% plots
from matplotlib.lines import Line2D
# only osmFISH_Ziesel (very similar to osmFISH_AllenSSp and osmFISH_AllenVISp
plt.style.use('ggplot')
plt.figure(figsize=(8, 2.5))
plt.boxplot([Starmap_before, Ziesel_before,MERFISH_before,seqFISH_before],
            positions = [1,4,7,10],boxprops=dict(color='red'),
            capprops=dict(color='red'),whiskerprops=dict(color='red'),flierprops=dict(color='red', markeredgecolor='red'),
            medianprops=dict(color='red'))
plt.boxplot([Starmap_after, Ziesel_after,MERFISH_after,seqFISH_after],
            positions = [2,5,8,11],boxprops=dict(color='blue'),
            capprops=dict(color='blue'),whiskerprops=dict(color='blue'),flierprops=dict(color='blue', markeredgecolor='blue'),
            medianprops=dict(color='blue'))
y = Starmap_before
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Starmap_after
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Ziesel_before
x = np.random.normal(4, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = Ziesel_after
x = np.random.normal(5, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = MERFISH_before
x = np.random.normal(7, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = MERFISH_after
x = np.random.normal(8, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = seqFISH_before
x = np.random.normal(10, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)
y = seqFISH_after
x = np.random.normal(11, 0.05, len(y))
plt.plot(x, y, 'k.', alpha=0.2)

plt.xticks((1.5,4.5,7.5,10.5),('STARmap_AllenVISp','osmFISH_Zeisel'
                          ,'MERFISH_Moffit','seqFISH_AllenVISp'),size=10)
plt.yticks(size=8)
plt.gca().set_ylim([-0.5,1.1])
plt.ylabel('Cosine similarity',size=12)
colors = ['red', 'blue']
lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
labels = ['Before PRECISE', 'After PRECISE']
plt.legend(lines, labels)
plt.show()
