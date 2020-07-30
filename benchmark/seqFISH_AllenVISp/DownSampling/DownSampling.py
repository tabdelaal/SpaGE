import os
os.chdir('seqFISH_AllenVISp/')

import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import scipy.stats as st
import sys
sys.path.insert(1,'Scripts/SpaGE/')
from principal_vectors import PVComputation

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

Gene_Order = np.intersect1d(seqFISH_data.columns,RNA_data.columns)

### SpaGE
SpaGE_imputed = pd.read_csv('Results/SpaGE_LeaveOneOut.csv',header=0,index_col=0,sep=',')
#SpaGE_New = pd.read_csv('Results/SpaGE_New_genes.csv',header=0,index_col=0,sep=',')

SpaGE_imputed = SpaGE_imputed.loc[:,Gene_Order]

SpaGE_seqCorr = pd.Series(index = Gene_Order)
for i in Gene_Order:
    SpaGE_seqCorr[i] = st.spearmanr(seqFISH_data[i],SpaGE_imputed[i])[0]
SpaGE_seqCorr[np.isnan(SpaGE_seqCorr)] = 0

SpaGE_seqCorr.sort_values(ascending=False,inplace=True)
test_set = SpaGE_seqCorr.index[0:50]

Common_data = RNA_data[Gene_Order]
Common_data = Common_data.drop(columns=test_set)

Corr = np.corrcoef(Common_data.T)
#for i in range(0,Corr.shape[0]):
#    Corr[i,i]=0
    
#plt.hist(np.abs(np.reshape(Corr,-1)),bins=np.arange(0,1.05,0.05))
#plt.show()   
# 0.7
removed_genes = []
for i in range(0,Corr.shape[0]):
    for j in range(i+1,Corr.shape[0]):
        if(np.abs(Corr[i,j]) > 0.7):
            Vi = np.var(Common_data.iloc[:,i])
            Vj = np.var(Common_data.iloc[:,j])
            if(Vi > Vj):
                removed_genes.append(Common_data.columns[j])
            else:
                removed_genes.append(Common_data.columns[i])
removed_genes= np.unique(removed_genes)

Common_data = Common_data.drop(columns=removed_genes)
Variance = np.var(Common_data)
Variance.sort_values(ascending=False,inplace=True)
Variance = Variance.append(pd.Series(0,index=removed_genes))

#### Novel Genes Expression Patterns ####
genes_to_impute = test_set
for i in [10,30,50,100,200,500,1000,2000,5000,7000,len(Variance)]:
    print(i)
    Imp_New_Genes = pd.DataFrame(np.zeros((seqFISH_data.shape[0],len(genes_to_impute))),columns=genes_to_impute)
    
    if(i>=50):
        n_factors = 50
        n_pv = 50
    else:
        n_factors = i
        n_pv = i
    
    dim_reduction = 'pca'
    dim_reduction_target = 'pca'
    
    pv_FISH_RNA = PVComputation(
            n_factors = n_factors,
            n_pv = n_pv,
            dim_reduction = dim_reduction,
            dim_reduction_target = dim_reduction_target
    )
    
    source_data = RNA_data_scaled[Variance.index[0:i]]
    target_data = seqFISH_data_scaled[Variance.index[0:i]]
    
    pv_FISH_RNA.fit(source_data,target_data)
    
    S = pv_FISH_RNA.source_components_.T
    
    Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
    S = S[:,0:Effective_n_pv]
    
    RNA_data_t = source_data.dot(S)
    FISH_exp_t = target_data.dot(S)
        
    nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',
                            metric = 'cosine').fit(RNA_data_t)
    distances, indices = nbrs.kneighbors(FISH_exp_t)
     
    for j in range(0,seqFISH_data.shape[0]):
        weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
        weights = weights/(len(weights)-1)
        Imp_New_Genes.iloc[j,:] = np.dot(weights,RNA_data[genes_to_impute].iloc[indices[j,:][distances[j,:] < 1]])
        
    Imp_New_Genes.to_csv('Results/' + str(i) +'SpaGE_New_genes.csv')
