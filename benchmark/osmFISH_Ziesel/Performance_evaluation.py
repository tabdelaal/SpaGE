import os
os.chdir('osmFISH_Ziesel/')

import numpy as np
import pandas as pd
import pickle
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.stats as st
from sklearn.neighbors import NearestNeighbors

plt.style.use('dark_background')
cmap = cm.get_cmap('viridis',20)
def plot_Imputed(df,destination,X,Y,New = False):
    for i in df.columns:
        fig, axs = plt.subplots()
        axs.axis("off")
        cmap = df[i]
        cmap[cmap > np.percentile(cmap,99)] = np.percentile(cmap,99)
        axs.scatter(X,Y,s=1,c=cmap)
        if(New):
             plt.savefig(destination +'/New_' + i + '.pdf')
        else:
            plt.savefig(destination +'/' + i + '.pdf')
        plt.close()

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data = datadict['osmFISH_data']
osmFISH_meta = datadict['osmFISH_meta']
del datadict

all_centroids  = osmFISH_meta[['X','Y']]

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
    
Moran_Is = Moran_I(osmFISH_data,all_centroids)

Gene_Order = osmFISH_data.columns

Moran_Is = Moran_Is[Gene_Order]

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')
SpaGE_New = pd.read_csv('Results/SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr[i] = st.spearmanr(osmFISH_data[i],SpaGE_imputed[i])[0]

plot_Imputed(osmFISH_data,'Figures/Original',osmFISH_meta['X'],
             osmFISH_meta['Y'])
plot_Imputed(SpaGE_imputed,'Figures/SpaGE_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'])
plot_Imputed(SpaGE_New,'Figures/SpaGE_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'],New = True)

### gimVI
gimVI_imputed = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')
gimVI_New = pd.read_csv('Results/gimVI_New_genes.csv',header=0,index_col=0,sep=',')

gimVI_imputed = gimVI_imputed.loc[:,[x.upper() for x in np.array(Gene_Order,dtype='str')]]

gimVI_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr[i] = st.spearmanr(osmFISH_data[i],gimVI_imputed[str(i).upper()])[0]
gimVI_Corr[np.isnan(gimVI_Corr)] = 0

plot_Imputed(gimVI_imputed,'Figures/gimVI_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'])
plot_Imputed(gimVI_New,'Figures/gimVI_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'],New = True)   

### Seurat
Seurat_imputed = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T
Seurat_New = pd.read_csv('Results/Seurat_New_genes.csv',header=0,index_col=0,sep=',').T

Seurat_imputed = Seurat_imputed.loc[:,Gene_Order]

Seurat_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr[i] = st.spearmanr(osmFISH_data[i],Seurat_imputed[i])[0]

plot_Imputed(Seurat_imputed,'Figures/Seurat_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'])
plot_Imputed(Seurat_New,'Figures/Seurat_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'],New = True)

### Liger
Liger_imputed = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T
Liger_New = pd.read_csv('Results/Liger_New_genes.csv',header=0,index_col=0,sep=',').T

Liger_imputed = Liger_imputed.loc[:,Gene_Order]

Liger_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr[i] = st.spearmanr(osmFISH_data[i],Liger_imputed[i])[0]
Liger_Corr[np.isnan(Liger_Corr)] = 0

plot_Imputed(Liger_imputed,'Figures/Liger_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'])
plot_Imputed(Liger_New,'Figures/Liger_Predicted',osmFISH_meta['X'],
             osmFISH_meta['Y'],New = True)

### Comparison plots
plt.style.use('ggplot')
fig, ax = plt.subplots(figsize=(3.7, 5.5))
ax.boxplot([SpaGE_Corr,Seurat_Corr, Liger_Corr,gimVI_Corr])

y = SpaGE_Corr
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = Seurat_Corr
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = Liger_Corr
x = np.random.normal(3, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)
y = gimVI_Corr
x = np.random.normal(4, 0.05, len(y))
plt.plot(x, y, 'g.', alpha=0.2)

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
    im = ax.scatter(X, Y, s=25, c=cmap)
    im.set_cmap('seismic')
    im.set_clim(vmin=0)
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
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.scatter(X, Y, s=25)
    Corr = st.spearmanr(X,Y)[0]
    plt.text(np.mean(plt.gca().get_xlim()),np.min(plt.gca().get_ylim()),'%1.3f'%Corr,color='black',size=9)
    plt.xticks(size=8)
    plt.yticks(size=8)
    plt.axis('scaled')
    plt.gca().set_ylim([-0.5,1])
    plt.gca().set_xlim([0,1])
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
