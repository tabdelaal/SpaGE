import os
os.chdir('osmFISH_Ziesel/')

import numpy as np
import pandas as pd
import scipy.stats as st
import pickle
import loompy

## Read RNA seq data
RNA_exp = pd.read_csv('data/Zeisel/expression_mRNA_17-Aug-2014.txt',header=None,sep='\t')

RNA_data = RNA_exp.iloc[11:19984,:]
RNA_data = RNA_data.drop([1],axis=1)
RNA_data.set_index(0,inplace=True)
RNA_data.columns = range(0,RNA_data.shape[1])
RNA_data = RNA_data.astype('float')

RNA_meta = RNA_exp.iloc[0:10,:]
RNA_meta = RNA_meta.drop([0],axis=1)
RNA_meta.set_index(1,inplace=True)
RNA_meta.columns = range(0,RNA_meta.shape[1])
RNA_meta = RNA_meta.T
del RNA_exp

RNA_data = RNA_data.loc[:,RNA_meta['tissue']=='sscortex']
RNA_meta = RNA_meta.loc[RNA_meta['tissue']=='sscortex',:]

RNA_data.columns = range(0,RNA_data.shape[1])
RNA_meta.index = range(0,RNA_meta.shape[0])

# filter lowely expressed genes
Genes_count = np.sum(RNA_data > 0, axis=1)
RNA_data = RNA_data.loc[Genes_count >=10,:]
del Genes_count

def Log_Norm_cpm(x):
    return np.log(((x/np.sum(x))*1000000) + 1)

RNA_data = RNA_data.apply(Log_Norm_cpm,axis=0)
RNA_data_scaled = pd.DataFrame(data=st.zscore(RNA_data.T),index = RNA_data.columns,columns=RNA_data.index)

datadict = dict()
datadict['RNA_data'] = RNA_data.T
datadict['RNA_data_scaled'] = RNA_data_scaled
datadict['RNA_meta'] = RNA_meta

with open('data/SpaGE_pkl/Ziesel.pkl','wb') as f:
    pickle.dump(datadict, f)

## Read loom data
ds = loompy.connect('data/osmFISH/osmFISH_SScortex_mouse_all_cells.loom')

FISH_Genes = ds.ra['Gene']
    
colAtr = ds.ca.keys()

df = pd.DataFrame()
for i in colAtr:
    df[i] = ds.ca[i]

osmFISH_meta = df.iloc[np.where(df.Valid == 1)[0], :]
osmFISH_data = ds[:,:]
osmFISH_data = osmFISH_data[:,np.where(df.Valid == 1)[0]]
osmFISH_data = pd.DataFrame(data= osmFISH_data, index= FISH_Genes)

del ds, colAtr, i, df, FISH_Genes

Cortex_Regions = ['Layer 2-3 lateral', 'Layer 2-3 medial', 'Layer 3-4', 
                  'Layer 4','Layer 5', 'Layer 6', 'Pia Layer 1']
Cortical = np.stack(i in Cortex_Regions for i in osmFISH_meta.Region)

osmFISH_meta = osmFISH_meta.iloc[Cortical,:]
osmFISH_data = osmFISH_data.iloc[:,Cortical]

cell_count = np.sum(osmFISH_data,axis=0)
def Log_Norm(x):
    return np.log(((x/np.sum(x))*np.median(cell_count)) + 1)

osmFISH_data = osmFISH_data.apply(Log_Norm,axis=0)
osmFISH_data_scaled = pd.DataFrame(data=st.zscore(osmFISH_data.T),index = osmFISH_data.columns,columns=osmFISH_data.index)

datadict = dict()
datadict['osmFISH_data'] = osmFISH_data.T
datadict['osmFISH_data_scaled'] = osmFISH_data_scaled
datadict['osmFISH_meta'] = osmFISH_meta

with open('data/SpaGE_pkl/osmFISH_Cortex.pkl','wb') as f:
    pickle.dump(datadict, f)
