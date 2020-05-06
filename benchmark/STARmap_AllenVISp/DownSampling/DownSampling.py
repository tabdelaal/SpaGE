import os
os.chdir('STARmap_AllenVISp/')

import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import sys
sys.path.insert(1,'Scripts/SpaGE/')
from principal_vectors import PVComputation

with open ('data/SpaGE_pkl/Starmap.pkl', 'rb') as f:
    datadict = pickle.load(f)

Starmap_data = datadict['Starmap_data']
Starmap_data_scaled = datadict['Starmap_data_scaled']
coords = datadict['coords']
del datadict

with open ('data/SpaGE_pkl/Allen_VISp.pkl', 'rb') as f:
    datadict = pickle.load(f)
    
RNA_data = datadict['RNA_data']
RNA_data_scaled = datadict['RNA_data_scaled']
del datadict

all_centroids  = np.vstack([c.mean(0) for c in coords])

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

Moran_Is.sort_values(ascending=False,inplace=True)
test_set = Moran_Is.index[0:50]

Common_data = RNA_data[np.intersect1d(Starmap_data.columns,RNA_data.columns)]
Common_data = Common_data.drop(columns=test_set)

Corr = np.corrcoef(Common_data.T)
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

genes_to_impute = test_set
for i in [10,30,50,100,200,len(Variance)]:
    Imp_New_Genes = pd.DataFrame(np.zeros((Starmap_data.shape[0],len(genes_to_impute))),columns=genes_to_impute)
    
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
    target_data = Starmap_data_scaled[Variance.index[0:i]]
    
    pv_FISH_RNA.fit(source_data,target_data)
    
    S = pv_FISH_RNA.source_components_.T
    
    Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
    S = S[:,0:Effective_n_pv]
    
    RNA_data_t = source_data.dot(S)
    FISH_exp_t = target_data.dot(S)
        
    nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',
                            metric = 'cosine').fit(RNA_data_t)
    distances, indices = nbrs.kneighbors(FISH_exp_t)
     
    for j in range(0,Starmap_data.shape[0]):
        weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
        weights = weights/(len(weights)-1)
        Imp_New_Genes.iloc[j,:] = np.dot(weights,RNA_data[genes_to_impute].iloc[indices[j,:][distances[j,:] < 1]])
        
    Imp_New_Genes.to_csv('Results/' + str(i) +'SpaGE_New_genes.csv')
