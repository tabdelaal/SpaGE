import os
os.chdir('seqFISH_AllenVISp/')


import numpy as np
import pandas as pd
import pickle
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
#from matplotlib import cm
import scipy.stats as st

with open ('data/SpaGE_pkl/seqFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

seqFISH_data = datadict['seqFISH_data']
seqFISH_meta= datadict['seqFISH_meta']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data = datadict['RNA_data']
del datadict

Gene_Order = np.intersect1d(seqFISH_data.columns,RNA_data.columns)

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_seqCorr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_seqCorr[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed[i])[0]
SpaGE_seqCorr[np.isnan(SpaGE_seqCorr)] = 0

os.chdir('STARmap_AllenVISp/')

with open ('data/SpaGE_pkl/Starmap.pkl', 'rb') as f:
    datadict = pickle.load(f)

coords = datadict['coords']
Starmap_data = datadict['Starmap_data']
del datadict

Gene_Order = np.intersect1d(Starmap_data.columns,RNA_data.columns)

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut_cutoff.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_starCorr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_starCorr[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed[i])[0]

def Compare_Correlations(X,Y):
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(X, Y, s=1)        
    ax.axvline(linestyle='--',color='gray')
    ax.axhline(linestyle='--',color='gray')
    plt.gca().set_ylim([-0.5,1])
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  
        np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k-')
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.xticks(size=8)
    plt.yticks(size=8)
    plt.show()

Starmap_seq_genes = np.intersect1d(Starmap_data.columns,seqFISH_data.columns)
Compare_Correlations(SpaGE_starCorr[Starmap_seq_genes],SpaGE_seqCorr[Starmap_seq_genes])
plt.xlabel('Spearman Correlation STARmap',size=12)
plt.ylabel('Spearman Correlation seqFISH',size=12)
plt.show()

fig, ax = plt.subplots(figsize=(3.7, 4.5))
ax.boxplot([SpaGE_starCorr[Starmap_seq_genes],SpaGE_seqCorr[Starmap_seq_genes]],widths=0.5)

y = SpaGE_starCorr[Starmap_seq_genes]
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)
y = SpaGE_seqCorr[Starmap_seq_genes]
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)

plt.xticks((1,2),('STARmap','seqFISH'),size=12)
plt.yticks(size=8)
plt.gca().set_ylim([-0.4,0.8])
plt.ylabel('Spearman Correlation',size=12)
#ax.set_aspect(aspect=3)
_,p_val = st.wilcoxon(SpaGE_starCorr[Starmap_seq_genes],SpaGE_seqCorr[Starmap_seq_genes])
plt.text(2,np.max(plt.gca().get_ylim()),'%1.2e'%p_val,color='black',size=8)
plt.show()

os.chdir('osmFISH_AllenVISp/')

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data = datadict['osmFISH_data']
del datadict

Gene_Order = osmFISH_data.columns

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut_cutoff.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_osmCorr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_osmCorr[i] = st.spearmanr(osmFISH_data[i],SpaGE_imputed[i])[0]

def Compare_Correlations(X,Y):
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(X, Y, s=25)        
    ax.axvline(linestyle='--',color='gray')
    ax.axhline(linestyle='--',color='gray')
    plt.gca().set_ylim([-0.5,1])
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  
        np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k-')
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.xticks(size=8)
    plt.yticks(size=8)
    plt.show()

osm_seq_genes = np.intersect1d(osmFISH_data.columns,seqFISH_data.columns)
Compare_Correlations(SpaGE_osmCorr[osm_seq_genes],SpaGE_seqCorr[osm_seq_genes])
plt.xlabel('Spearman Correlation osmFISH',size=12)
plt.ylabel('Spearman Correlation seqFISH',size=12)
plt.show()

fig, ax = plt.subplots(figsize=(3.7, 4.5))
ax.boxplot([SpaGE_osmCorr[osm_seq_genes],SpaGE_seqCorr[osm_seq_genes]],widths=0.5)

y = SpaGE_osmCorr[osm_seq_genes]
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = SpaGE_seqCorr[osm_seq_genes]
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)

plt.xticks((1,2),('osmFISH','seqFISH'),size=12)
plt.yticks(size=8)
plt.gca().set_ylim([-0.5,1])
plt.ylabel('Spearman Correlation',size=12)
#ax.set_aspect(aspect=3)
_,p_val = st.wilcoxon(SpaGE_osmCorr[osm_seq_genes],SpaGE_seqCorr[osm_seq_genes])
plt.text(2,np.max(plt.gca().get_ylim()),'%1.2e'%p_val,color='black',size=8)
plt.show()
