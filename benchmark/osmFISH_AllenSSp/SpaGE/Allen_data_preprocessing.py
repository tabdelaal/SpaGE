import os
os.chdir('osmFISH_AllenSSp/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle

RNA_data = pd.read_csv('data/Allen_SSp/SSp_exons_matrix.csv',
                       header=0,index_col=0,sep=',')

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

with open('data/SpaGE_pkl/Allen_SSp.pkl','wb') as f:
    pickle.dump(datadict, f)
