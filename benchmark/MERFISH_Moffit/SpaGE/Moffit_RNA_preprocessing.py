import os
os.chdir('MERFISH_Moffit/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle
import scipy.io as io

genes = pd.read_csv('data/Moffit_RNA/GSE113576/genes.tsv',sep='\t',header=None)
barcodes = pd.read_csv('data/Moffit_RNA/GSE113576/barcodes.tsv',sep='\t',header=None)

genes = np.array(genes.loc[:,1])
barcodes = np.array(barcodes.loc[:,0])
RNA_data = io.mmread('data/Moffit_RNA/GSE113576/matrix.mtx')
RNA_data = RNA_data.todense()
RNA_data = pd.DataFrame(RNA_data,index=genes,columns=barcodes)

# filter lowely expressed genes
Genes_count = np.sum(RNA_data > 0, axis=1)
RNA_data = RNA_data.loc[Genes_count >=10,:]
del Genes_count

def Log_Norm(x):
    return np.log(((x/np.sum(x))*1000000) + 1)

RNA_data = RNA_data.apply(Log_Norm,axis=0)
RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.T),index = RNA_data.columns,columns=RNA_data.index)

datadict = dict()
datadict['RNA_data'] = RNA_data.T
datadict['RNA_data_scaled'] = RNA_data_scaled

with open('data/SpaGE_pkl/Moffit_RNA.pkl','wb') as f:
    pickle.dump(datadict, f, protocol=4)