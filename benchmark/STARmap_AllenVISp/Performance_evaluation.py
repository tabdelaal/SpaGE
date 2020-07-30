import os
os.chdir('STARmap_AllenVISp/')

import numpy as np
import pandas as pd
import pickle
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import scipy.stats as st
from sklearn.neighbors import NearestNeighbors
from matplotlib import cm

with open ('data/SpaGE_pkl/Starmap.pkl', 'rb') as f:
    datadict = pickle.load(f)

coords = datadict['coords']
Starmap_data = datadict['Starmap_data']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data = datadict['RNA_data']
del datadict

all_centroids  = np.vstack([c.mean(0) for c in coords])

plt.style.use('dark_background')
cmap = cm.get_cmap('viridis',20)

def Moran_I(SpatialData,XYmap):
    XYnbrs = NearestNeighbors(n_neighbors=5, algorithm='auto',metric = 'euclidean').fit(XYmap)
    XYdistances, XYindices = XYnbrs.kneighbors(XYmap)
    W = np.zeros((SpatialData.shape[0],SpatialData.shape[0]))
    for i in range(0,SpatialData.shape[0]):
        W[i,XYindices[i,:]]=1

    for i in range(0,SpatialData.shape[0]):
        W[i,i]=0
    
    I = pd.Series(index=SpatialData.columns)
    for k in SpatialData.columns:
        X_minus_mean = np.array(SpatialData[k] - np.mean(SpatialData[k]))
        X_minus_mean = np.reshape(X_minus_mean,(len(X_minus_mean),1))
        Nom = np.sum(np.multiply(W,np.matmul(X_minus_mean,X_minus_mean.T)))
        Den = np.sum(np.multiply(X_minus_mean,X_minus_mean))
        I[k] = (len(SpatialData[k])/np.sum(W))*(Nom/Den)
    return(I)
    
Moran_Is = Moran_I(Starmap_data,all_centroids)
    
Gene_Order = np.intersect1d(Starmap_data.columns,RNA_data.columns)

Moran_Is = Moran_Is[Gene_Order]

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr[i] = st.spearmanr(Starmap_data[i],SpaGE_imputed[i])[0]
    
### gimVI
gimVI_imputed = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')

gimVI_imputed.columns = Gene_Order

gimVI_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr[i] = st.spearmanr(Starmap_data[i],gimVI_imputed[i])[0]
gimVI_Corr[np.isnan(gimVI_Corr)] = 0

### Seurat
Seurat_imputed = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Seurat_imputed = Seurat_imputed.loc[:,Gene_Order]
Seurat_imputed.index = range(0,Seurat_imputed.shape[0])
cell_labels = pd.read_csv('data/Starmap/visual_1020/20180505_BY3_1kgenes/class_labels.csv',
                          header=0,sep=',')
Starmap_data_Seurat = Starmap_data.drop(np.where(cell_labels['ClusterName']=='HPC')[0],axis=0)

Seurat_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr[i] = st.spearmanr(Starmap_data_Seurat[i],Seurat_imputed[i])[0]

### Liger
Liger_imputed = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Liger_imputed = Liger_imputed.loc[:,Gene_Order]
Liger_imputed.index = range(0,Liger_imputed.shape[0])

Liger_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr[i] = st.spearmanr(Starmap_data[i],Liger_imputed[i])[0]
Liger_Corr[np.isnan(Liger_Corr)] = 0

### Comparison plots
plt.style.use('ggplot')
fig, ax = plt.subplots(figsize=(3.7, 5.5))
ax.boxplot([SpaGE_Corr,Seurat_Corr, Liger_Corr,gimVI_Corr])

y = SpaGE_Corr
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)
y = Seurat_Corr
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)
y = Liger_Corr
x = np.random.normal(3, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)
y = gimVI_Corr
x = np.random.normal(4, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=1, alpha=0.2)

plt.xticks((1,2,3,4),('SpaGE', 'Seurat', 'Liger','gimVI'),size=12)
plt.yticks(size=8)
plt.gca().set_ylim([-0.5,1])
plt.ylabel('Spearman Correlation',size=12)
ax.set_aspect(aspect=3)
_,p_val = st.wilcoxon(SpaGE_Corr,Seurat_Corr)
plt.text(2,np.max(plt.gca().get_ylim()),'%1.2e'%p_val,color='black',size=8)
_,p_val = st.wilcoxon(SpaGE_Corr,Liger_Corr)
plt.text(3,np.max(plt.gca().get_ylim()),'%1.2e'%p_val,color='black',size=8)
_,p_val = st.wilcoxon(SpaGE_Corr,gimVI_Corr)
plt.text(4,np.max(plt.gca().get_ylim()),'%1.2e'%p_val,color='black',size=8)
plt.show()

def Compare_Correlations(X,Y):
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    cmap = Moran_Is
    ax.axvline(linestyle='--',color='gray')
    ax.axhline(linestyle='--',color='gray')
    im = ax.scatter(X, Y, s=1, c=cmap)
    im.set_cmap('seismic')
    plt.gca().set_ylim([-0.5,1])
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  
        np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, 'k-')
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.xticks((-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),size=8)
    plt.yticks((-0.4,-0.2,0,0.2,0.4,0.6,0.8,1),size=8)
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=8) 
    cbar.ax.set_ylabel("Moran's I statistic",fontsize=12)
    plt.show()

Compare_Correlations(Seurat_Corr,SpaGE_Corr)
plt.xlabel('Spearman Correlation Seurat',size=12)
plt.ylabel('Spearman Correlation SpaGE',size=12)
plt.show()

Compare_Correlations(Liger_Corr,SpaGE_Corr)
plt.xlabel('Spearman Correlation Liger',size=12)
plt.ylabel('Spearman Correlation SpaGE',size=12)
plt.show()

Compare_Correlations(gimVI_Corr,SpaGE_Corr)
plt.xlabel('Spearman Correlation gimVI',size=12)
plt.ylabel('Spearman Correlation SpaGE',size=12)
plt.show()

def Correlation_vs_Moran(X,Y):
    fig, ax = plt.subplots(figsize=(4.8, 4.8))
    ax.scatter(X, Y, s=1)
    Corr = st.spearmanr(X,Y)[0]
    plt.text(np.mean(plt.gca().get_xlim()),np.min(plt.gca().get_ylim()),'%1.3f'%Corr,color='black',size=9)
    plt.xticks(size=8)
    plt.yticks(size=8)
    plt.axis('scaled')
    plt.gca().set_ylim([-0.5,1])
    plt.gca().set_xlim([-0.2,1])
    plt.show()


Correlation_vs_Moran(Moran_Is,SpaGE_Corr)
plt.xlabel("Moran's I",size=12)
plt.ylabel('Spearman Correlation SpaGE',size=12)
plt.show()

Correlation_vs_Moran(Moran_Is,Seurat_Corr)
plt.xlabel("Moran's I",size=12)
plt.ylabel('Spearman Correlation Seurat',size=12)
plt.show()

Correlation_vs_Moran(Moran_Is,Liger_Corr)
plt.xlabel("Moran's I",size=12)
plt.ylabel('Spearman Correlation Liger',size=12)
plt.show()

Correlation_vs_Moran(Moran_Is,gimVI_Corr)
plt.xlabel("Moran's I",size=12)
plt.ylabel('Spearman Correlation gimVI',size=12)
plt.show()
