import os
os.chdir('osmFISH_AllenVISp/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

RNA_data = pd.read_csv('data/Allen_VISp/mouse_VISp_2018-06-14_exon-matrix.csv',
                       header=0,index_col=0,sep=',')

Genes = pd.read_csv('data/Allen_VISp/mouse_VISp_2018-06-14_genes-rows.csv',
                        header=0,sep=',')
RNA_data.index = Genes.gene_symbol
del Genes

# filter lowely expressed genes
Genes_count = np.sum(RNA_data > 0, axis=1)
RNA_data = RNA_data.loc[Genes_count >=10,:]
del Genes_count

# filter low quality cells
meta_data = pd.read_csv('data/Allen_VISp/mouse_VISp_2018-06-14_samples-columns.csv',
                        header=0,sep=',')
HighQualityCells = (meta_data['class'] != 'No Class') & (meta_data['class'] != 'Low Quality')
RNA_data = RNA_data.iloc[:,np.where(HighQualityCells)[0]]
del meta_data, HighQualityCells

def Log_Norm(x):
    return np.log(((x/np.sum(x))*1000000) + 1)

RNA_data = RNA_data.apply(Log_Norm,axis=0)
RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.T),index = RNA_data.columns,columns=RNA_data.index)

datadict = dict()
datadict['RNA_data'] = RNA_data.T
datadict['RNA_data_scaled'] = RNA_data_scaled

with open('data/SpaGE_pkl/Allen_VISp.pkl','wb') as f:
    pickle.dump(datadict, f)