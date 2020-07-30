import os
os.chdir('MERFISH_Moffit/')

import numpy as np
import pandas as pd
import pickle
import matplotlib
matplotlib.use('qt5agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import scipy.stats as st

with open ('data/SpaGE_pkl/MERFISH.pkl', 'rb') as f:
    datadict = pickle.load(f)

MERFISH_data = datadict['MERFISH_data']
del datadict

with open ('data/SpaGE_pkl/Moffit_RNA.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data = datadict['RNA_data']
del datadict

Gene_Order = np.intersect1d(MERFISH_data.columns,RNA_data.columns)

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_Corr[i] = st.spearmanr(MERFISH_data[i],SpaGE_imputed[i])[0]
    
### gimVI
gimVI_imputed = pd.read_csv('Results/gimVI_LeaveOneOut.csv',header=0,index_col=0,sep=',')
gimVI_imputed = gimVI_imputed.drop(columns='AVPR2')

gimVI_imputed = gimVI_imputed.loc[:,[x.upper() for x in np.array(Gene_Order,dtype='str')]]

gimVI_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    gimVI_Corr[i] = st.spearmanr(MERFISH_data[i],gimVI_imputed[str(i).upper()])[0]
gimVI_Corr[np.isnan(gimVI_Corr)] = 0


### Seurat
Seurat_imputed = pd.read_csv('Results/Seurat_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Seurat_imputed = Seurat_imputed.loc[:,Gene_Order]

Seurat_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Seurat_Corr[i] = st.spearmanr(MERFISH_data[i],Seurat_imputed[i])[0]

### Liger
Liger_imputed = pd.read_csv('Results/Liger_LeaveOneOut.csv',header=0,index_col=0,sep=',').T

Liger_imputed = Liger_imputed.loc[:,Gene_Order]

Liger_Corr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    Liger_Corr[i] = st.spearmanr(MERFISH_data[i],Liger_imputed[i])[0]
Liger_Corr[np.isnan(Liger_Corr)] = 0

### Comparison plots
plt.style.use('ggplot')
fig, ax = plt.subplots(figsize=(3.7, 5.5))
ax.boxplot([SpaGE_Corr,Seurat_Corr, Liger_Corr,gimVI_Corr])

y = SpaGE_Corr
x = np.random.normal(1, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=3, alpha=0.2)
y = Seurat_Corr
x = np.random.normal(2, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=3, alpha=0.2)
y = Liger_Corr
x = np.random.normal(3, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=3, alpha=0.2)
y = gimVI_Corr
x = np.random.normal(4, 0.05, len(y))
plt.plot(x, y, 'g.', markersize=3, alpha=0.2)

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
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.scatter(X, Y, s=5)        
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
