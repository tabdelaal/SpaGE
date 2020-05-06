import os
os.chdir('osmFISH_Ziesel/')

import pickle
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
import time as tm

with open ('data/SpaGE_pkl/Ziesel.pkl', 'rb') as f:
    datadict = pickle.load(f)

RNA_data = datadict['RNA_data']
RNA_data_scaled = datadict['RNA_data_scaled']
del datadict

with open ('data/SpaGE_pkl/osmFISH_Cortex.pkl', 'rb') as f:
    datadict = pickle.load(f)

osmFISH_data = datadict['osmFISH_data']
osmFISH_data_scaled = datadict['osmFISH_data_scaled']
osmFISH_meta= datadict['osmFISH_meta']
del datadict

#### Leave One Out Validation ####
Common_data = RNA_data_scaled[np.intersect1d(osmFISH_data_scaled.columns,RNA_data_scaled.columns)]
Imp_Genes = pd.DataFrame(columns=Common_data.columns)
precise_time = []
knn_time = []
for i in Common_data.columns:
    print(i)
    start = tm.time()
    from principal_vectors import PVComputation

    n_factors = 30
    n_pv = 30
    dim_reduction = 'pca'
    dim_reduction_target = 'pca'

    pv_FISH_RNA = PVComputation(
            n_factors = n_factors,
            n_pv = n_pv,
            dim_reduction = dim_reduction,
            dim_reduction_target = dim_reduction_target
    )

    pv_FISH_RNA.fit(Common_data.drop(i,axis=1),osmFISH_data_scaled[Common_data.columns].drop(i,axis=1))

    S = pv_FISH_RNA.source_components_.T
    
    Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
    S = S[:,0:Effective_n_pv]

    Common_data_t = Common_data.drop(i,axis=1).dot(S)
    FISH_exp_t = osmFISH_data_scaled[Common_data.columns].drop(i,axis=1).dot(S)
    precise_time.append(tm.time()-start)
        
    start = tm.time()
    nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',metric = 'cosine').fit(Common_data_t)
    distances, indices = nbrs.kneighbors(FISH_exp_t)
    
    Imp_Gene = np.zeros(osmFISH_data.shape[0])
    for j in range(0,osmFISH_data.shape[0]):
        weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
        weights = weights/(len(weights)-1)
        Imp_Gene[j] = np.sum(np.multiply(RNA_data[i][indices[j,:][distances[j,:] < 1]],weights))
    Imp_Gene[np.isnan(Imp_Gene)] = 0
    Imp_Genes[i] = Imp_Gene
    knn_time.append(tm.time()-start)

Imp_Genes.to_csv('Results/SpaGE_LeaveOneOut.csv')
precise_time = pd.DataFrame(precise_time)
knn_time = pd.DataFrame(knn_time)
precise_time.to_csv('Results/SpaGE_PreciseTime.csv', index = False)
knn_time.to_csv('Results/SpaGE_knnTime.csv', index = False)

#### Novel Genes Expression Patterns ####
Common_data = RNA_data_scaled[np.intersect1d(osmFISH_data_scaled.columns,RNA_data_scaled.columns)]
genes_to_impute = ["Tesc","Pvrl3","Grm2"]
Imp_New_Genes = pd.DataFrame(np.zeros((osmFISH_data.shape[0],len(genes_to_impute))),columns=genes_to_impute)

from principal_vectors import PVComputation

n_factors = 30
n_pv = 30
dim_reduction = 'pca'
dim_reduction_target = 'pca'

pv_FISH_RNA = PVComputation(
        n_factors = n_factors,
        n_pv = n_pv,
        dim_reduction = dim_reduction,
        dim_reduction_target = dim_reduction_target
)

pv_FISH_RNA.fit(Common_data,osmFISH_data_scaled[Common_data.columns])

S = pv_FISH_RNA.source_components_.T
    
Effective_n_pv = sum(np.diag(pv_FISH_RNA.cosine_similarity_matrix_) > 0.3)
S = S[:,0:Effective_n_pv]

Common_data_t = Common_data.dot(S)
FISH_exp_t = osmFISH_data_scaled[Common_data.columns].dot(S)
    
nbrs = NearestNeighbors(n_neighbors=50, algorithm='auto',
                        metric = 'cosine').fit(Common_data_t)
distances, indices = nbrs.kneighbors(FISH_exp_t)

for j in range(0,osmFISH_data.shape[0]):

    weights = 1-(distances[j,:][distances[j,:]<1])/(np.sum(distances[j,:][distances[j,:]<1]))
    weights = weights/(len(weights)-1)
    Imp_New_Genes.iloc[j,:] = np.dot(weights,RNA_data[genes_to_impute].iloc[indices[j,:][distances[j,:] < 1]])

Imp_New_Genes.to_csv('Results/SpaGE_New_genes.csv')